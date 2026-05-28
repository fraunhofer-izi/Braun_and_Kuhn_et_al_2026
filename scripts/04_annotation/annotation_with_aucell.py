import sys
import argparse
from pathlib import Path
import spatialdata as sd
import scanpy as sc
import logging
import json
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from spatial_transcriptomics_analysis.scripts.anndata_utils import annotate_with_custom_aucell_threshold,check_l2_matches_l1
from spatial_transcriptomics_analysis.scripts.anndata_plotting import crop_region_by_fraction

logger = logging.getLogger(__name__)
def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

# ------------ Argument Parsing ----------------
def parse_args():
    parser = argparse.ArgumentParser(description="Annoation of table based on gene sets based " \
    "of AUCell scores from processed ZARR files.")
    parser.add_argument(
        "--zarr-path",
        type=str,
        help="Path to .zarr file to process (overrides value in config if given)."
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="YAML config file with 'sample_name', and gene set parameters."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for figures and filtered zarrs."
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity level"
    )
    return parser.parse_args()

def map_to_l1_broad(adata, mapping_path: str):
    """Map adata.obs['L2_fits_L1'] (L2 subtypes) to broad L1 classes using a JSON mapping."""
    if mapping_path is None:
        logger.warning("No L2->L1 broad mapping_path provided; skipping L2_fits_L1_broad mapping.")
        return adata

    mapping_path = Path(mapping_path)
    if not mapping_path.exists():
        raise FileNotFoundError(f"L2->L1 broad mapping JSON not found: {mapping_path}")

    with open(mapping_path, "r") as f:
        l2_to_l1_raw = json.load(f)

    # reverse mapping: subcelltype -> broad cell type
    l2_to_l1 = {
        subcelltype: l1
        for l1, subcelltypes in l2_to_l1_raw.items()
        for subcelltype in subcelltypes
    }

    adata.obs["L2_fits_L1_broad"] = adata.obs["L2_fits_L1"].map(l2_to_l1)

    mask_unmapped = adata.obs["L2_fits_L1"].notna() & adata.obs["L2_fits_L1_broad"].isna()
    unmapped = adata.obs.loc[mask_unmapped, "L2_fits_L1"].unique()
    if len(unmapped) > 0:
        logger.warning("Unmapped L2 cell types found (no broad mapping): %s", unmapped)

    return adata

# ------------ Create Output Directories ----------------
def create_output_dirs(base_dir: Path):
    data_anno_dir = base_dir / "zarr_after_anno"
    figures_anno_dir = base_dir / "figures" / "annotation"
    data_anno_dir.mkdir(parents=True, exist_ok=True)
    figures_anno_dir.mkdir(parents=True, exist_ok=True)
    return data_anno_dir, figures_anno_dir

def main():
    args = parse_args()
    setup_logging(args.log_level)

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    zarr_path = Path(args.zarr_path).resolve() if args.zarr_path else Path(config["zarr_path"]).resolve()
    sample_name = config.get("sample_name", zarr_path.stem)
    output_dir = Path(args.output_dir).resolve()

    # Load QC thresholds from config or defaults
    LEVEL1_GENEMARKERS = config.get("l1_json_file")
    LEVEL1_CUTOFF = config.get("l1_auc_cutoff", 0.05)
    LEVEL2_GENEMARKERS = config.get("l2_json_file")
    LEVEL2_CUTOFF = config.get("l2_auc_cutoff", 0.1)
    LEVEL1_TO_LEVEL2_MAPPING = config.get("l2_to_l1_mapping_file")

    data_anno_dir, figures_anno_dir = create_output_dirs(output_dir)

    logger.info(f"\n🚀 Processing sample: {sample_name}")
    logger.info(f"📂 Loading from: {zarr_path}")

    sdata = sd.read_zarr(zarr_path)

    adata = sdata.tables["cell_bins"].copy()
    adata = annotate_with_custom_aucell_threshold(
        adata,
        json_file=LEVEL1_GENEMARKERS,
        score_prefix="custom_markers_l1",
        score_threshold=LEVEL1_CUTOFF
    )
    adata = annotate_with_custom_aucell_threshold(
        adata,
        json_file=LEVEL2_GENEMARKERS,
        score_prefix="custom_markers_l2",
        score_threshold=LEVEL2_CUTOFF
    )
    adata = check_l2_matches_l1(
        adata,
        json_file=LEVEL1_TO_LEVEL2_MAPPING,
        l2_col="custom_markers_l2_celltype_filtered",
        l1_passing_col="custom_markers_l1_celltypes_passing",
        output_col="L2_fits_L1",
        l1_from_l2_col="mapped_L2_to_L1",
        fits_boolean_col="L2_L1_consistent"
    )
    valid = adata.obs["L2_fits_L1"].notna().sum()
    invalid = adata.obs["L2_fits_L1"].isna().sum()

    adata = map_to_l1_broad(adata, mapping_path=LEVEL1_TO_LEVEL2_MAPPING)

    logger.info("Number of L2 celltypes that fit L1 mapping:", valid)
    logger.info("Number of mismatches (NaN):", invalid)
    
    sdata.tables["cell_bins"] = adata

    CROPPING_REGION = config.get("cropping_fraction", None)
    if CROPPING_REGION is not None:
        Y_FRAC = CROPPING_REGION.get("yfrac", None)
        X_FRAC = CROPPING_REGION.get("xfrac", None)
    crop_box = crop_region_by_fraction(adata, yfrac=Y_FRAC, xfrac=X_FRAC)

    sc.pl.spatial(
        adata,
        color="custom_markers_l1_celltype_filtered",
        crop_coord=crop_box,
        size=0.7,
        legend_fontsize=12,
        title=None,
        frameon=False,
    )
    plt.savefig(figures_anno_dir / f"{sample_name}_annotation_l1.png",
            dpi=300, bbox_inches="tight")
    plt.close()
    sc.pl.spatial(
        adata,
        color="custom_markers_l2_celltype_filtered",
        crop_coord=crop_box,
        size=0.7,
        legend_fontsize=12,
        title=None,
        frameon=False,
    )
    plt.savefig(figures_anno_dir / f"{sample_name}_annotation_l2.png",
            dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.spatial(
        adata,
        color="L2_fits_L1",
        crop_coord=crop_box,
        size=0.7,
        legend_fontsize=12,
        title=None,
        frameon=False,
    )
    plt.savefig(figures_anno_dir / f"{sample_name}_annotation_l1_fits_l2.png",
            dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.spatial(
        adata,
        color="L2_fits_L1_broad",
        crop_coord=crop_box,
        size=0.7,
        legend_fontsize=12,
        title=None,
        frameon=False,
    )
    plt.savefig(figures_anno_dir / f"{sample_name}_annotation_l2_mapped_to_l1.png",
            dpi=300, bbox_inches="tight")
    plt.close()

    sdata.write(data_anno_dir / f"{sample_name}.zarr", overwrite=True)
    logger.info("\n✅ Sample Annotation complete!")

if __name__ == "__main__":
    main()