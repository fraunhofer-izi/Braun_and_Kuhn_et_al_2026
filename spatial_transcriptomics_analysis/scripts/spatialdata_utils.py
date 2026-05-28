from typing import Dict
import numpy as np
import spatialdata as sd
from spatialdata.models import TableModel
import spatialdata_io
from spatialdata_io._constants._constants import VisiumHDKeys
import geopandas as gpd
import scanpy as sc
from geopandas import GeoDataFrame
import json
from pathlib import Path
from PIL import Image
from typing import Dict, Tuple

import logging
logger = logging.getLogger(__name__)

def read_visium_hd_segmented(outs_path: str, sample_id: str) -> TableModel:
    """
    Build a TableModel (AnnData) for Visium HD segmented outputs:
    - reads raw_feature_cell_matrix.h5
    - attaches spatial images + scalefactors under adata.uns['spatial']
    - sets region / instance keys for SpatialData integration
    - adds qc metrics + centroid-based spatial coordinates
    - estimates spot_diameter_fullres from mean cell area
    """
    outs_path = Path(outs_path)
    bin_size = "segmented_outputs"
    path_bin = outs_path / bin_size
    path_bin_spatial = path_bin / VisiumHDKeys.SPATIAL

    # Load gene expression
    counts_file = "raw_feature_cell_matrix.h5"
    adata = sc.read_10x_h5(path_bin / counts_file, gex_only=False)

    # Load scalefactors and images
    with open(path_bin_spatial / VisiumHDKeys.SCALEFACTORS_FILE) as f:
        scalefactors = json.load(f)

    hires_img = np.array(Image.open(path_bin_spatial / "tissue_hires_image.png"))
    lowres_img = np.array(Image.open(path_bin_spatial / "tissue_lowres_image.png"))

    library_id = "visium_hd_segmentation"
    adata.uns["spatial"] = {
        library_id: {
            "images": {
                "hires": hires_img,
                "lowres": lowres_img,
            },
            "scalefactors": scalefactors,
            "metadata": {
                "source_image_path": "tissue_hires_image.png"
            }
        }
    }

    # Assign region info
    shapes_name = "cell_segmentation"
    adata.obs[VisiumHDKeys.INSTANCE_KEY] = np.arange(len(adata))
    adata.obs[VisiumHDKeys.REGION_KEY] = shapes_name
    adata.obs[VisiumHDKeys.REGION_KEY] = adata.obs[VisiumHDKeys.REGION_KEY].astype("category")

    # Attach QC metrics (optional)
    adata.var_names_make_unique()
    adata = adata[:, ~adata.var_names.str.startswith("DEPRECATED")].copy()
    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito", "ribo"], inplace=True)

    # Load cell shapes
    shapes = spatialdata_io.geojson(outs_path / "segmented_outputs" / "cell_segmentations.geojson", coordinate_system=sample_id)

    # Match shape order to adata
    centroids = shapes.geometry.centroid
    adata.obsm["spatial"] = np.vstack([centroids.x.values, centroids.y.values]).T

    # Estimate spot diameter
    areas = shapes.geometry.area
    mean_area = np.mean(areas)
    diameter_pixels = 2 * np.sqrt(mean_area / np.pi)
    adata.uns["spatial"][library_id]["scalefactors"]["spot_diameter_fullres"] = diameter_pixels
    
    return TableModel.parse(
        adata=adata,
        region=shapes_name,
        region_key=str(VisiumHDKeys.REGION_KEY),
        instance_key=str(VisiumHDKeys.INSTANCE_KEY),
    )


def load_visiumhd_scalefactors_and_images(outs_path: Path) -> Tuple[dict, Dict[str, np.ndarray], str]:
    """
    Load scalefactors and required images from VisiumHD output folder.
    """
    spatial_dir = outs_path / "segmented_outputs" / VisiumHDKeys.SPATIAL

    with open(spatial_dir / VisiumHDKeys.SCALEFACTORS_FILE) as f:
        scalefactors = json.load(f)

    images = {
        "hires": np.array(Image.open(spatial_dir / "tissue_hires_image.png")),
        "lowres": np.array(Image.open(spatial_dir / "tissue_lowres_image.png")),
    }

    return scalefactors, images, str(spatial_dir / "tissue_hires_image.png")


def add_table_to_shapes(
    adata: sc.AnnData,
    shapes: GeoDataFrame,
    shapes_name: str,
    library_id: str,
    scalefactors: dict,
    images: Dict[str, np.ndarray],
    source_image_path: str
) -> TableModel:
    """
    Link AnnData and Segmentation shapes, add spatial metadata for scanpy/squidpy plotting and return a TableModel.
    """
    adata.uns["spatial"] = {
        library_id: {
            "images": images,
            "scalefactors": scalefactors,
            "metadata": {"source_image_path": source_image_path}
        }
    }

    adata.obs["location_id"] = np.arange(len(adata))
    adata.obs["region"] = shapes_name
    adata.obs["region"] = adata.obs["region"].astype("category")

    centroids = shapes.geometry.centroid
    adata.obsm["spatial"] = np.vstack([centroids.x.values, centroids.y.values]).T

    areas = shapes.geometry.area
    diameter_pixels = 2 * np.sqrt(np.mean(areas) / np.pi)
    adata.uns["spatial"][library_id]["scalefactors"]["spot_diameter_fullres"] = diameter_pixels

    return TableModel.parse(
        adata=adata,
        region=shapes_name,
        region_key="region",
        instance_key="location_id",
    )

def prepare_sample_adatas(sample_paths):
    """
    Reads SpatialData once and extracts 'cell_bins' AnnData per sample.
    Returns a dict of sample_name -> AnnData and sample_name -> SpatialData.
    """
    adatas = {}
    sdatas = {}
    for sample_name, path in sample_paths.items():
        sdata = sd.read_zarr(path)
        if "cell_bins" not in sdata.tables:
            raise KeyError(f"'cell_bins' not found in {sample_name}")
        adata = sdata.tables["cell_bins"].copy()
        sc.pp.calculate_qc_metrics(adata, inplace=True)  # ensures all relevant metrics are present
        adatas[sample_name] = adata
        sdatas[sample_name] = sdata
    return adatas, sdatas