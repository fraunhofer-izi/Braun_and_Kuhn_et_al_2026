# Copyright 2026 Fraunhofer-Gesellschaft zur Förderung der angewandten
# Forschung e.V.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

# Logger konfigurieren
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def calculate_qc(adata):
    """
    Calculate quality control (QC) metrics for an AnnData object.

    Notes:
    ------
    - Mitochondrial genes are identified by names starting with "MT-".
    - Ribosomal genes are identified by names starting with "RPS" or "RPL".
    - The `sc.pp.calculate_qc_metrics` function from Scanpy is used to compute
      the QC metrics.
    """
    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito", "ribo"], inplace=True)
    return adata


def filter_adata_genes_cells(adata, min_cells, min_counts):
    """Filters genes and cells in an AnnData object based on minimum cell and count thresholds.

    Args:
        adata (anndata.AnnData): The annotated data matrix to be filtered.
        min_cells (int): Minimum number of cells a gene must be expressed in to be retained.
        min_counts (int): Minimum number of counts a cell must have to be retained.

    Return:
        anndata.AnnData: The filtered AnnData object with genes and cells removed according to the specified thresholds.
    """
    logging.info("Filter AnnData based on min_cells and min_counts")
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    return adata


def is_outlier(adata, metric: str, nmads: int):
    """
    Identify outliers in a given metric of an AnnData object based on the median and
    median absolute deviation (MAD).

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix. The metric to evaluate should be present in `adata.obs`.
    metric : str
        The key in `adata.obs` corresponding to the metric to evaluate for outliers.
    nmads : int
        The number of median absolute deviations (MADs) to use as the threshold for
        identifying outliers.

    Returns:
    --------
    pandas.Series
        A boolean series indicating whether each observation is an outlier (True) or not (False).
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def qc_plots(adata, output_path):
    """
    Generate quality control (QC) plots for single-cell RNA sequencing data.

    This function creates a 2x2 grid of plots to visualize various quality control metrics
    from the AnnData object and saves the resulting figure to the specified output path.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix. The function expects the following columns in `adata.obs`:
        - "total_counts": Total counts per cell.
        - "pct_counts_mito": Percentage of mitochondrial gene counts per cell.
        - "n_genes_by_counts": Number of genes detected per cell.
    output_path : str
        File path where the generated QC plot will be saved.

    Plots:
    ------
    1. Histogram of total counts.
    2. Violin plot for mitochondrial content percentage (`pct_counts_mito`).
    3. Scatter plot of total counts vs. number of genes by counts, colored by mitochondrial content.
    4. Violin plot for the number of genes by counts (`n_genes_by_counts`).

    The resulting figure is saved as a high-resolution image (300 dpi) at the specified output path.

    Returns:
    --------
    None
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid

    # 1. Histogram of total counts
    sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0, 0])
    axes[0, 0].set_title("Total Counts Histogram")
    axes[0, 0].set_xlabel("total_counts")

    # 2. Violin plot for pct_counts_mito
    sns.violinplot(
        y=adata.obs["pct_counts_mito"], ax=axes[0, 1], inner="box", color="skyblue"
    )
    axes[0, 1].set_title("Mitochondrial Content (%)")
    axes[0, 1].set_ylabel("pct_counts_mito")

    # 3. Scatter plot with color bar
    sc = axes[1, 0].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mito"],
        cmap="viridis",
        s=2,
    )
    cbar = fig.colorbar(sc, ax=axes[1, 0])
    cbar.set_label("pct_counts_mito")
    axes[1, 0].set_title("Total Counts vs N Genes by Counts")
    axes[1, 0].set_xlabel("total_counts")
    axes[1, 0].set_ylabel("n_genes_by_counts")

    # 4. Violin for n_genes_by_counts
    sns.violinplot(
        y=adata.obs["n_genes_by_counts"], ax=axes[1, 1], inner="box", color="lightgreen"
    )
    axes[1, 1].set_title("Number of Genes by Counts")
    axes[1, 1].set_ylabel("n_genes_by_counts")

    # Final layout
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)  # You can specify any path here
    plt.show()


def qc_plots_combined(adata, filtered_adata, output_path):
    """
    Compare quality control metrics before and after filtering,
    including a bar plot of total cell numbers.
    """
    fig, axes = plt.subplots(3, 2, figsize=(14, 16))

    # Combine metadata into one DataFrame with a label
    # Combine metadata into one DataFrame with a label
    # Combine metadata safely with source labels, dropping the index to avoid reindex issues
    adata_df = adata.obs[
        ["total_counts", "pct_counts_mito", "n_genes_by_counts"]
    ].copy()
    adata_df["source"] = "before"
    filtered_df = filtered_adata.obs[
        ["total_counts", "pct_counts_mito", "n_genes_by_counts"]
    ].copy()
    filtered_df["source"] = "after"

    # Drop index explicitly to ensure no duplicate labels
    combined = pd.concat([adata_df, filtered_df], ignore_index=True)

    # 1. Histogram of total counts
    sns.histplot(
        data=combined,
        x="total_counts",
        hue="source",
        bins=100,
        ax=axes[0, 0],
        element="step",
        stat="density",
        common_norm=False,
    )
    axes[0, 0].set_title("Total Counts Histogram")
    axes[0, 0].set_xlabel("total_counts")

    # 2. Violin plot for pct_counts_mito
    sns.violinplot(
        data=combined,
        y="pct_counts_mito",
        x="source",
        ax=axes[0, 1],
        inner="box",
        palette="pastel",
    )
    axes[0, 1].set_title("Mitochondrial Content (%)")
    axes[0, 1].set_ylabel("pct_counts_mito")

    # 3. Scatter plot with color bar
    sc = axes[1, 0].scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mito"],
        cmap="viridis",
        s=2,
        label="before",
    )
    axes[1, 0].scatter(
        filtered_adata.obs["total_counts"],
        filtered_adata.obs["n_genes_by_counts"],
        c=filtered_adata.obs["pct_counts_mito"],
        cmap="cool",
        s=2,
        label="after",
        alpha=0.6,
    )
    cbar = fig.colorbar(sc, ax=axes[1, 0])
    cbar.set_label("pct_counts_mito")
    axes[1, 0].set_title("Total Counts vs N Genes by Counts")
    axes[1, 0].set_xlabel("total_counts")
    axes[1, 0].set_ylabel("n_genes_by_counts")
    axes[1, 0].legend()

    # 4. Violin plot for n_genes_by_counts
    sns.violinplot(
        data=combined,
        y="n_genes_by_counts",
        x="source",
        ax=axes[1, 1],
        inner="box",
        palette="Set2",
    )
    axes[1, 1].set_title("Number of Genes by Counts")
    axes[1, 1].set_ylabel("n_genes_by_counts")

    # 5. Bar plot: Number of cells before vs. after
    cell_counts = pd.DataFrame(
        {"source": ["before", "after"], "n_cells": [adata.n_obs, filtered_adata.n_obs]}
    )
    sns.barplot(
        data=cell_counts, x="source", y="n_cells", ax=axes[2, 0], palette="muted"
    )
    axes[2, 0].set_title("Total Number of Cells")
    axes[2, 0].set_ylabel("Cell Count")
    axes[2, 0].set_xlabel("")

    # Hide the unused subplot (bottom right)
    axes[2, 1].axis("off")

    # Final layout
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
