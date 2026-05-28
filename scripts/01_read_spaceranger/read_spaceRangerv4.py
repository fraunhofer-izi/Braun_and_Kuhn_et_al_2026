from pathlib import Path
import scanpy as sc
import numpy as np
from PIL import Image
import spatialdata_io
import argparse
import logging
import yaml

from spatial_transcriptomics_analysis.scripts.spatialdata_plotting import (
    crop_sdata_to_shape_bounds
)
from spatial_transcriptomics_analysis.scripts.spatialdata_utils import (
    read_visium_hd_segmented
)
logger = logging.getLogger(__name__)
def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

def parse_args():
    parser = argparse.ArgumentParser(description="Read VisiumHD sample with cell segmentation.")
    parser.add_argument(
        "--spaceranger-outs",
        type=str,
        help="Path to the output directory of a SpaceRanger run."
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="YAML config file with 'sample_name'"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Base output directory"
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity level"
    )
    return parser.parse_args()

def create_output_dirs(base_dir: Path):
    zarr_dir = Path(base_dir) / "zarr_before_qc"
    zarr_dir.mkdir(parents=True, exist_ok=True)
    return zarr_dir


def main():
    args = parse_args()
    setup_logging(args.log_level)
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    outs_path = Path(args.spaceranger_outs).resolve() if args.spaceranger_outs else Path(config["sample_path"]).resolve()
    sample_id = config.get("sample_name")
    zarr_dir = create_output_dirs(args.output_dir)

    logger.info("💾 Read VisiumHD output")
    # Create base SpatialData object
    sdata = spatialdata_io.visium_hd(outs_path, sample_id)

    # Add segmentation shapes
    logger.info("💾 Add Segmentation")
    sdata["nucleus_segmentation"] = spatialdata_io.geojson(
        outs_path / "segmented_outputs" / "nucleus_segmentations.geojson",
        coordinate_system=sample_id
    )
    sdata["cell_segmentation"] = spatialdata_io.geojson(
        outs_path / "segmented_outputs" / "cell_segmentations.geojson",
        coordinate_system=sample_id
    )
    # Add table
    logger.info("💾 Add AnnData")
    sdata["cell_bins"] = read_visium_hd_segmented(outs_path, sample_id)
    sdata_crop = crop_sdata_to_shape_bounds(sdata, f"{sample_id}_square_016um", sample_id)
    logger.info("💾 Saving Zarr output")

    output_path = zarr_dir / f"{sample_id}.zarr"
    sdata_crop.write(output_path, overwrite=True)
    logger.info(f"💾 Saved Zarr output to {output_path}")

if __name__ == "__main__":
    main()
