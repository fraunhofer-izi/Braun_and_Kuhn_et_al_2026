import argparse
from pathlib import Path
import logging
import yaml
import spatialdata as sd
import scanpy as sc
import spatialdata_plot
import matplotlib.pyplot as plt
import time
import numpy as np

from spatial_transcriptomics_analysis.scripts.anndata_utils import (
    cluster_based_on_gex,
    clustering_with_hergast
)

logger = logging.getLogger(__name__)
def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

# ---------------- Argument Parsing ----------------

def parse_args():
    parser = argparse.ArgumentParser(description="Run clustering on a sample .zarr")
    parser.add_argument(
        "--zarr-path",
        type=str,
        help="Path to .zarr file to process (overrides value in config if given)."
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

# ---------------- Directory Setup ----------------

def create_output_dirs(base_dir: Path):
    clusters_dir = base_dir / "zarr_after_cluster"
    figs_dir = base_dir / "figures" / "clusters"
    clusters_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)
    return clusters_dir, figs_dir

# ---------------- Main Logic ----------------

def main():
    args = parse_args()
    setup_logging(args.log_level)
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    zarr_path = Path(args.zarr_path).resolve() if args.zarr_path else Path(config["zarr_path"]).resolve()
    sample_name = config.get("sample_name", zarr_path.stem)
    output_dir = Path(args.output_dir).resolve()
    clusters_dir, figs_dir = create_output_dirs(output_dir)

    logger.info(f"🚀 Clustering sample: {sample_name}")
    sdata = sd.read_zarr(zarr_path)

    for name, table in sdata.tables.items():
        if name in ["square_002um", "square_016um", "square_008um"]:
            continue
        logger.info(f"🧪 Clustering table: {name}")

        adata = table.copy()
        adata.raw = adata.copy() #save raw counts
        sc.settings.n_jobs = -1
        adata = cluster_based_on_gex(adata)
        adata = clustering_with_hergast(adata)
        sdata.tables[name] = adata

    # Save clustered ZARR
    sdata.write(clusters_dir / f"{sample_name}.zarr", overwrite=True)
    logger.info("✅ Clustering done & saved!")

if __name__ == "__main__":
    main()
