import sys
import argparse
from pathlib import Path
import spatialdata as sd
import spatialdata_plot
import scanpy as sc
import logging
import pandas as pd
import yaml

from spatial_transcriptomics_analysis.scripts.anndata_qc import calculate_qc, qc_plots_combined

logger = logging.getLogger(__name__)
def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )
# ------------ Argument Parsing ----------------
def parse_args():
    parser = argparse.ArgumentParser(description="Run QC and visualization from processed ZARR files.")
    parser.add_argument(
        "--zarr-path",
        type=str,
        help="Path to .zarr file to process (overrides value in config if given)."
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="YAML config file with 'sample_name', and optional QC parameters."
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

# ------------ Create Output Directories ----------------
def create_output_dirs(base_dir: Path):
    data_qc_dir = base_dir / "zarr_after_qc"
    figures_qc_dir = base_dir / "figures" / "qc"
    data_qc_dir.mkdir(parents=True, exist_ok=True)
    figures_qc_dir.mkdir(parents=True, exist_ok=True)
    return data_qc_dir, figures_qc_dir

def main():
    args = parse_args()
    setup_logging(args.log_level)

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    zarr_path = Path(args.zarr_path).resolve() if args.zarr_path else Path(config["zarr_path"]).resolve()
    sample_name = config.get("sample_name", zarr_path.stem)
    output_dir = Path(args.output_dir).resolve()

    # Load QC thresholds from config or defaults
    MIN_GENES_08 = config.get("min_genes_08", 30)
    MIN_COUNTS_08 = config.get("min_counts_08", 50)
    MIN_GENES = config.get("min_genes", 100)
    MIN_COUNTS = config.get("min_counts", 200)
    MAX_MITO = config.get("max_mito", 15)

    data_qc_dir, figures_qc_dir = create_output_dirs(output_dir)

    logger.info(f"\n🚀 Processing sample: {sample_name}")
    logger.info(f"📂 Loading from: {zarr_path}")

    sdata = sd.read_zarr(zarr_path)

    for name, table in sdata.tables.items():
        qc_numbers = {"before_qc": {}, "after_qc": {}}
        if name in ["square_002um"]:
            continue

        logger.info(f"🧪 QC on table: {name}")
        adata = table.copy()
        adata.obs_names_make_unique()

        calculate_qc(adata)
        qc_numbers["before_qc"]["n_cells"] = adata.n_obs

        shape_name = f"{sample_name}_{name}" if "square" in name else "cell_segmentation"
        sdata.tables[name] = adata
        sdata.pl.render_images(f"{sample_name}_hires_image")\
                 .pl.render_shapes(shape_name, color="total_counts", cmap="viridis", table_name=name)\
                 .pl.show(coordinate_systems=sample_name, dpi=600,
                          save=figures_qc_dir / f"{sample_name}_{name}_ncounts_before_QC.png")

        filtered = adata[adata.obs["pct_counts_mito"] < MAX_MITO].copy()
        filtered.obs_names_make_unique()

        if name == "square_016um":
            sc.pp.filter_cells(filtered, min_counts=MIN_COUNTS)
            sc.pp.filter_cells(filtered, min_genes=MIN_GENES)
        else:
            sc.pp.filter_cells(filtered, min_counts=MIN_COUNTS_08)
            sc.pp.filter_cells(filtered, min_genes=MIN_GENES_08)
        sc.pp.filter_genes(filtered, min_cells=5)

        sdata.tables[name] = filtered
        sdata.pl.render_images(f"{sample_name}_hires_image")\
                 .pl.render_shapes(shape_name, color="total_counts", cmap="viridis", table_name=name)\
                 .pl.show(coordinate_systems=sample_name, dpi=600,
                          save=figures_qc_dir / f"{sample_name}_{name}_ncounts_after_gene_cell_removal_QC.png")

        # if name == "cell_bins":
        #     cell_shapes = sdata.shapes["cell_segmentation"]
        #     microns_per_pixel = filtered.uns["spatial"]["visium_hd_segmentation"]["scalefactors"]["microns_per_pixel"]
        #     cell_shapes["area_um2"] = cell_shapes.geometry.area * (microns_per_pixel ** 2)

        #     filtered_shapes = cell_shapes[(cell_shapes["area_um2"] > 10) & (cell_shapes["area_um2"] < 400)]
        #     valid_ids = filtered_shapes.index
        #     tmp = filtered[filtered.obs["location_id"].isin(valid_ids)].copy()

        #     sdata.shapes["cell_segmentation"] = filtered_shapes
        #     filtered = tmp
        
        sdata.tables[name] = filtered
        qc_numbers["after_qc"]["n_cells"] = filtered.n_obs
        df = pd.DataFrame.from_dict(qc_numbers, orient="index")
        df.index.name = "stage"
        df.to_csv(figures_qc_dir / f"{sample_name}_{name}_qc_numbers.csv")

        qc_plots_combined(adata, filtered, figures_qc_dir / f"{sample_name}_{name}_after_qc.png")

        sdata.pl.render_images(f"{sample_name}_hires_image")\
                 .pl.render_shapes(shape_name, color="total_counts", cmap="viridis", table_name=name)\
                 .pl.show(coordinate_systems=sample_name, dpi=600,
                          save=figures_qc_dir / f"{sample_name}_{name}_ncounts_after_cellsizeremoval_QC.png")

    sdata.write(data_qc_dir / f"{sample_name}.zarr", overwrite=True)
    logger.info("\n✅ Sample QC complete!")

if __name__ == "__main__":
    main()