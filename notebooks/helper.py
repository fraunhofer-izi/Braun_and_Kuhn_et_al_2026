# Copyright 2026 Fraunhofer-Gesellschaft zur Förderung der angewandten
# Forschung e.V.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License.

import json
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib.colors as mcolors
import liana as li
import seaborn as sns
import math
from liana.method import (
    cellchat,
    cellphonedb,
    connectome,
    geometric_mean,
    logfc,
    natmi,
    singlecellsignalr,
)
from liana.mt import rank_aggregate
from matplotlib.patches import Patch
from plotnine import (
    element_blank,
    element_text,
    theme,
)
from spatial_transcriptomics_analysis.scripts.anndata_plotting import (
    crop_region_by_fraction,
)


def load_configs(config_dir):
    """
    Load all YAML config files from a directory.

    Parameters
    ----------
    config_dir : str or Path
        Directory containing YAML config files.

    Returns
    -------
    dict
        Dictionary with sample_name as keys and config dicts as values.
    """

    config_dir = Path(config_dir)

    configs = {}

    for config_path in config_dir.glob("*.yaml"):

        with open(config_path, "r") as f:
            cfg = yaml.safe_load(f)

        sample = cfg["sample_name"]

        configs[sample] = cfg

    return configs


def get_configs(sample_name, configs):
    """
    Extract plotting-related config values for a sample.
    """

    cfg = configs[sample_name]

    crop = cfg.get("cropping_fraction", None)
    cluster_key = cfg.get("clustering_group")
    cluster_palette = cfg.get("cluster_color_mapping", None)
    color_mappings = cfg.get("color_mappings")

    return crop, cluster_key, cluster_palette, color_mappings


# ---------------------------------------------------------
# Add CAR_pos column
# ---------------------------------------------------------
def add_car_positive_column(adata, car_genes=["Cilta", "ciltacel"]):
    """
    Adds a boolean 'CAR_pos' column based on expression of known CAR transgenes
    """

    # Case-insensitive matching in var_names
    matches = [
        g for g in adata.var_names if g.lower() in [c.lower() for c in car_genes]
    ]

    if matches:
        car_gene = matches[0]

        Xg = adata[:, car_gene].X
        if hasattr(Xg, "toarray"):
            expr = Xg.toarray().ravel()
        else:
            expr = np.array(Xg).ravel()

        adata.obs["CAR_pos"] = expr > 0
    else:
        adata.obs["CAR_pos"] = False
        print("No CAR transgene found in var_names (Cilta / ciltacel).")


# ---------------------------------------------------------
# Add T cell column
# ---------------------------------------------------------
def assign_tcell(row, T_CELL_TYPES):
    ct = row["L2_fits_L1"]
    return ct in T_CELL_TYPES  # True / False


# ---------------------------------------------------------
# Assign CAR status for T cells only
# ---------------------------------------------------------
def assign_tcell_car(row):
    if row["T_cell"]:
        return "CAR_pos" if row["CAR_pos"] else "CAR_neg"
    else:
        return "non_T_cell"


# ---------------------------------------------------------
# Assign CAR status within T-cell subtype
# ---------------------------------------------------------
def assign_tsubtype_car(row, T_CELL_TYPES):
    subtype = row["L2_fits_L1"]
    if subtype in T_CELL_TYPES:
        return f"{subtype}_CAR_pos" if row["CAR_pos"] else f"{subtype}_CAR_neg"
    else:
        return "non_T_cell"


def save_cluster_legend(
    l1_colors,
    outfile,
    exclude_labels=None,
    figsize=(2.6, 0.9),
    fontsize=8,
    ncol=1,
    dpi=300,
):
    """
    Create and save a standalone legend figure.

    Parameters
    ----------
    l1_colors : dict
        Dictionary mapping labels to colors.
    outfile : str
        Path to save the legend image.
    exclude_labels : list, optional
        Labels to exclude from the legend.
    figsize : tuple, optional
        Figure size.
    fontsize : int, optional
        Font size for legend text.
    ncol : int, optional
        Number of legend columns.
    dpi : int, optional
        DPI for saved figure.
    """

    if exclude_labels is None:
        exclude_labels = ["Undetermined", "Ambiguous", "T and NK cells"]

    # Filter labels
    filtered_colors = {k: v for k, v in l1_colors.items() if k not in exclude_labels}

    labels = list(filtered_colors.keys())

    handles = [
        Patch(facecolor=filtered_colors[label], edgecolor="none", label=label)
        for label in labels
    ]

    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.axis("off")

    fig.legend(
        handles=handles,
        labels=labels,
        loc="center",
        ncol=ncol,
        frameon=False,
        handlelength=1.0,
        handleheight=0.8,
        handletextpad=0.4,
        columnspacing=1.0,
        fontsize=fontsize,
    )

    fig.savefig(
        outfile,
        dpi=dpi,
        bbox_inches="tight",
        transparent=True,
    )

    plt.close(fig)


def plot_spatial_regions(
    sample,
    roi_presents,
    roi_name,
    adatas,
    configs,
    cmap,
    figdir_base="./results/publication_plots",
    cluster_size=0.7,
):
    """
    Generate publication-ready spatial plots for a sample ROI.

    Parameters
    ----------
    sample : str
        Sample name.

    roi_presents: dict
        Dictionary containing ROI presets for samples.

    roi_name : str
        Name of ROI defined in roi_presents.

    adatas : dict
        Dictionary of AnnData objects.

    configs : dict
        Config dictionary loaded with load_configs().

    cmap : matplotlib colormap
        Colormap used for plotting.

    figdir_base : str
        Base directory for output plots.

    cluster_size : float
        Spot size for spatial plots.
    """

    # --------------------------------------------------------------
    # Validate ROI
    # --------------------------------------------------------------
    if sample not in roi_presents:
        raise ValueError(f"No ROI presets found for sample '{sample}'")

    if roi_name not in roi_presents[sample]:
        raise ValueError(f"ROI '{roi_name}' not found for sample '{sample}'")

    # --------------------------------------------------------------
    # Load ROI coordinates
    # --------------------------------------------------------------
    roi = roi_presents[sample][roi_name]

    xfrac = roi["x"]
    yfrac = roi["y"]

    # --------------------------------------------------------------
    # Output directory
    # --------------------------------------------------------------
    figdir = os.path.join(
        figdir_base,
        sample,
        "spatial_regions",
        roi_name,
    )

    os.makedirs(figdir, exist_ok=True)

    sc.settings.figdir = figdir

    # --------------------------------------------------------------
    # Load sample configs
    # --------------------------------------------------------------
    crop, cluster_key, cluster_palette, color_mappings = get_configs(
        sample,
        configs,
    )

    # --------------------------------------------------------------
    # Load color mappings
    # --------------------------------------------------------------
    if isinstance(color_mappings, str):

        with open(color_mappings, "r") as f:
            color_mappings = json.load(f)

    l1_colors = color_mappings["cell_type_colors"].copy()

    # --------------------------------------------------------------
    # Add CAR colors
    # --------------------------------------------------------------
    l1_colors.update(
        {
            "CAR_pos": "#70FFF1",
            "CAR_neg": "#3E6EFF",
            "NK": "#CC3363",
        }
    )

    # --------------------------------------------------------------
    # Load adata
    # --------------------------------------------------------------
    adata = adatas[sample].copy()

    # --------------------------------------------------------------
    # Generate crop box
    # --------------------------------------------------------------
    crop_box = crop_region_by_fraction(
        adata,
        yfrac=yfrac,
        xfrac=xfrac,
    )

    # --------------------------------------------------------------
    # Full spatial clusters
    # --------------------------------------------------------------
    sc.pl.spatial(
        adata,
        color=cluster_key,
        size=cluster_size,
        palette=cluster_palette,
        cmap=cmap,
        frameon=False,
        show=True,
        save=f"_{roi_name}_clusters_full.png",
    )

    # --------------------------------------------------------------
    # Cropped spatial clusters
    # --------------------------------------------------------------
    sc.pl.spatial(
        adata,
        color=cluster_key,
        size=cluster_size,
        crop_coord=crop_box,
        palette=cluster_palette,
        cmap=cmap,
        frameon=False,
        show=True,
        save=f"_{roi_name}_clusters_crop.png",
    )

    # --------------------------------------------------------------
    # CAR status plot
    # --------------------------------------------------------------
    sc.pl.spatial(
        adata,
        color="L2_fits_L1_and_CAR_status",
        size=cluster_size,
        crop_coord=crop_box,
        cmap=cmap,
        palette=l1_colors,
        frameon=False,
        show=True,
        save=f"_{roi_name}_CAR_status_crop.png",
    )

    # --------------------------------------------------------------
    # Save annotation legend
    # --------------------------------------------------------------
    legend_outfile = os.path.join(
        figdir,
        f"{sample}_{roi_name}_annotation_legend.png",
    )

    save_cluster_legend(
        l1_colors=l1_colors,
        outfile=legend_outfile,
    )

    # --------------------------------------------------------------
    # Save cluster legend
    # --------------------------------------------------------------
    legend_outfile = os.path.join(
        figdir,
        f"{sample}_{roi_name}_cluster_legend.png",
    )

    save_cluster_legend(
        l1_colors=cluster_palette,
        outfile=legend_outfile,
        ncol=len(cluster_palette),
    )

    print(f"Finished plotting sample '{sample}' " f"for ROI '{roi_name}'")


def plot_celltype_composition(
    sample,
    adata,
    color_map_file,
    celltype_column="L2_fits_L1_broad",
    tcell_column="Tsubtype_CAR_status",
    non_tcell_label="non_T_cell",
    figsize=(4 / 2.54, 3 / 2.54),
    wspace=1.0,
    fontsize=6,
    figdir_base="./results/publication_plots",
    save_name="celltype_composition.png",
):
    """
    Plot stacked barplots for:
    1. All cell types
    2. T-cell phenotypes

    Parameters
    ----------
    sample : str
        Sample name.
    adata : AnnData
        AnnData object for the sample.
    color_map_file : str
        Path to color_maps.json.
    celltype_column : str
        Column for broad cell types.
    tcell_column : str
        Column for T-cell phenotypes.
    non_tcell_label : str
        Label used for non-T cells.
    figsize : tuple
        Figure size.
    wspace : float
        Space between subplots.
    fontsize : int
        Title fontsize.
    """

    # ----------------------------------------------------------
    # Load color maps
    # ----------------------------------------------------------
    with open(color_map_file) as f:
        color_maps = json.load(f)

    cell_type_colors = color_maps["cell_type_colors"]
    t_cell_colors = color_maps.get("t_cell_colors", {})

    # ----------------------------------------------------------
    # Plot 1: all cell types
    # ----------------------------------------------------------
    counts1 = adata.obs[celltype_column].value_counts()

    df1 = pd.DataFrame([counts1])

    colors1 = [cell_type_colors.get(cat, "lightgrey") for cat in df1.columns]

    # ----------------------------------------------------------
    # Plot 2: T-cell phenotypes
    # ----------------------------------------------------------
    adata_T = adata[adata.obs[tcell_column] != non_tcell_label]

    counts2 = adata_T.obs[tcell_column].value_counts()

    df2 = pd.DataFrame([counts2])

    colors2 = [t_cell_colors.get(cat, "lightgrey") for cat in df2.columns]

    # ----------------------------------------------------------
    # Figure
    # ----------------------------------------------------------
    fig, axes = plt.subplots(
        1,
        2,
        figsize=figsize,
        gridspec_kw={"wspace": wspace},
    )

    # ----------------------------------------------------------
    # Plot 1
    # ----------------------------------------------------------
    df1.plot(
        kind="bar",
        stacked=True,
        ax=axes[0],
        width=0.8,
        color=colors1,
        legend=False,
        edgecolor="none",
    )

    axes[0].set_title(
        "All celltypes",
        fontsize=fontsize,
    )

    axes[0].set_xticks([])
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Counts", fontsize=6)
    axes[0].tick_params(axis="y", labelsize=5)

    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    # ----------------------------------------------------------
    # Plot 2
    # ----------------------------------------------------------
    df2.plot(
        kind="bar",
        stacked=True,
        ax=axes[1],
        width=0.8,
        color=colors2,
        legend=False,
        edgecolor="none",
    )

    axes[1].set_title(
        "(CAR) T cell\nphenotypes",
        fontsize=fontsize,
    )

    axes[1].set_xticks([])
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    axes[1].tick_params(axis="y", labelsize=5)

    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    # ----------------------------------------------------------
    # Grid styling
    # ----------------------------------------------------------
    for ax in axes:

        ax.grid(
            axis="y",
            color="lightgrey",
            linewidth=0.3,
            alpha=0.4,
        )

    plt.tight_layout()
    sample_dir = os.path.join(figdir_base, sample)
    os.makedirs(sample_dir, exist_ok=True)

    plt.savefig(
        os.path.join(sample_dir, save_name),
        bbox_inches="tight",
        dpi=300,
    )
    plt.show()

    return


def plot_top_genes_per_cluster(
    adata,
    sample,
    number_genes=5,
    column="leiden_HERGAST_0_4",
    figdir_base="./results/publication_plots",
):

    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    top_genes = []

    for group in groups:
        names = result["names"][group]
        lfc = result["logfoldchanges"][group]
        padj = result["pvals_adj"][group]

        # keep only significant positive markers
        all_genes = [
            (g, fc, p) for g, fc, p in zip(names, lfc, padj) if fc > 0 and p < 0.05
        ]

        # select top genes
        selected = all_genes[:number_genes]

        print(f"Cluster {group}: " f"Selected top genes: {[g[0] for g in selected]}")

        top_genes.extend([g[0] for g in selected])

    # remove duplicates while preserving order
    unique_genes = list(dict.fromkeys(top_genes))

    # ---------------------------------------------------------
    # SAVE LOGIC
    # ---------------------------------------------------------

    figdir = os.path.join(
        figdir_base,
        sample,
        "top_genes_per_cluster",
    )

    os.makedirs(figdir, exist_ok=True)

    sc.settings.figdir = figdir

    # ---------------------------------------------------------
    # PLOTS
    # ---------------------------------------------------------

    sc.pl.dotplot(
        adata,
        var_names=unique_genes,
        groupby=column,
        standard_scale="var",
        save=f"_{sample}dotplot.png",
        show=True,
    )

    sc.pl.matrixplot(
        adata,
        var_names=unique_genes,
        groupby=column,
        standard_scale="var",
        save=f"_{sample}matrixplot.png",
        show=False,
    )


def prepare_tcell_composition_data(
    adata,
    cluster_key,
    plot_mode="fraction",
    cluster_label_map=None,
):
    """
    Prepare T-cell composition dataframe for plotting.
    """

    # ----------------------------------------------------------
    # Total cells per region
    # ----------------------------------------------------------
    total_cells = adata.obs[cluster_key].astype(str).value_counts().sort_index()

    # ----------------------------------------------------------
    # T-cell subtype counts
    # ----------------------------------------------------------
    t_counts = pd.crosstab(
        adata.obs[cluster_key].astype(str),
        adata.obs["Tsubtype_CAR_status"],
    )

    t_counts = t_counts.drop(
        columns=[
            "none",
            "unknown",
            "unassigned",
            "non_T_cell",
        ],
        errors="ignore",
    )

    # ----------------------------------------------------------
    # Normalize
    # ----------------------------------------------------------
    norm = t_counts.div(total_cells, axis=0)

    # ----------------------------------------------------------
    # Plot mode
    # ----------------------------------------------------------
    if plot_mode == "fraction":

        plot_df = norm.copy()

        xlabel = "Fraction of all cells"
        title = "T cell composition per region"
        outfile_suffix = "fractions"

    elif plot_mode == "counts":

        plot_df = t_counts.copy()

        xlabel = "Absolute T-cell counts"
        title = "Absolute T-cell counts per region"
        outfile_suffix = "counts"

    else:
        raise ValueError("plot_mode must be 'fraction' or 'counts'")

    # ----------------------------------------------------------
    # Relabel clusters
    # ----------------------------------------------------------
    if cluster_label_map is not None:

        plot_df.index = [
            cluster_label_map.get(int(i), i) if str(i).isdigit() else i
            for i in plot_df.index
        ]

        desired_order = list(cluster_label_map.values())

        desired_order = [lab for lab in desired_order if lab in plot_df.index]

        plot_df = plot_df.reindex(desired_order)

    return plot_df, xlabel, title, outfile_suffix


def plot_tcell_composition(
    sample_id,
    adata,
    plot_df,
    color_map_file,
    comparison_dir,
    outfile_suffix,
    xlabel,
    title,
    figsize=(6 / 2.54, 4 / 2.54),
    xlabel_fontsize=8,
    ylabel_fontsize=8,
    legend_fontsize=7,
):
    """
    Plot T-cell composition stacked horizontal barplot.
    """

    # ----------------------------------------------------------
    # Load colors
    # ----------------------------------------------------------
    with open(color_map_file) as f:
        color_maps = json.load(f)

    color_map = color_maps.get("t_cell_colors", {})

    categories = plot_df.columns.tolist()

    colors = [color_map.get(c, "#CCCCCC") for c in categories]

    # ----------------------------------------------------------
    # Figure
    # ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)

    # ----------------------------------------------------------
    # Plot
    # ----------------------------------------------------------
    plot_df.plot(
        kind="barh",
        stacked=True,
        color=colors,
        ax=ax,
        legend=False,
    )

    # ----------------------------------------------------------
    # Legend
    # ----------------------------------------------------------
    ax.legend(
        categories,
        title="T cell subtype",
        loc="upper center",
        bbox_to_anchor=(1.6, 1.2),
        ncol=1,
        frameon=False,
        fontsize=legend_fontsize,
        title_fontsize=legend_fontsize,
    )

    # ----------------------------------------------------------
    # Styling
    # ----------------------------------------------------------
    ax.set_xlabel(
        xlabel,
        fontsize=xlabel_fontsize,
    )

    ax.set_title(
        title,
        fontsize=xlabel_fontsize,
    )

    ax.set_ylabel("")

    ax.tick_params(
        axis="x",
        labelsize=xlabel_fontsize,
    )

    ax.tick_params(
        axis="y",
        labelsize=ylabel_fontsize,
    )

    ax.invert_yaxis()

    ax.grid(
        True,
        axis="x",
        linewidth=0.5,
        alpha=0.3,
    )

    ax.grid(
        True,
        axis="y",
        linewidth=0.5,
        alpha=0.3,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    # ----------------------------------------------------------
    # Save
    # ----------------------------------------------------------
    figdir = os.path.join(comparison_dir, sample_id)

    os.makedirs(figdir, exist_ok=True)

    sc.settings.figdir = figdir

    outfile = os.path.join(
        figdir,
        f"Tcell_composition_{outfile_suffix}.png",
    )

    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()


def neighbor_stats_per_cell(adata, cluster_key, source, targets):
    """
    For a given source cell type, compute:
    1) average number of neighbors per cell (absolute)
    2) neighbor composition (fraction / percentage)

    Returns a DataFrame with both metrics.
    """
    A = adata.obsp["spatial_connectivities"]  # Graph (CSR sparse)
    labels = adata.obs[cluster_key].astype("category")
    cats = labels.cat.categories
    codes = labels.cat.codes.to_numpy()  # 0..n_cls-1

    src = int(np.where(cats == source)[0][0])
    src_idx = np.where(codes == src)[0]

    if len(src_idx) == 0:
        raise ValueError(f"No cells found for source: {source}")

    results = []

    total_edges = 0
    edge_counts = {}

    # first pass: count edges
    for target in targets:
        tgt = int(np.where(cats == target)[0][0])
        tgt_idx = np.where(codes == tgt)[0]

        sub = A[src_idx][:, tgt_idx]
        n_edges = sub.sum()

        edge_counts[target] = n_edges
        total_edges += n_edges

    # second pass: compute per-cell and composition
    for target in targets:
        n_edges = edge_counts[target]

        neighbors_per_cell = n_edges / len(src_idx)
        composition = n_edges / total_edges if total_edges > 0 else np.nan

        results.append(
            {
                "source": source,
                "target": target,
                "neighbors_per_cell": neighbors_per_cell,
                "composition_fraction": composition,
                "composition_percent": composition * 100,
            }
        )

    return pd.DataFrame(results)


def autopct_absolute(values, fmt="{:.1f}"):
    """
    Display absolute numbers (e.g., neighbors_per_cell) inside the pie.
    The pie slice sizes are still determined by the provided percentages.
    """
    values = list(values)

    def _autopct(pct):
        total = sum(values)
        val = pct * total / 100.0
        return fmt.format(val)

    return _autopct


def get_colors_for_targets(targets, palette, fallback="#CCCCCC"):
    """Map a list of labels to colors using a palette dict."""
    return [palette.get(t, fallback) for t in targets]


def annotate_pie(ax, wedges, values, base_radius=1.35, small_thresh=3.1):
    """
    Annotate pie wedges with values outside the pie.
    Small values are pushed further out and get leader lines to avoid overlap.
    """
    values = list(values)

    for w, val in zip(wedges, values):
        ang = (w.theta2 + w.theta1) / 2.0
        ang_rad = np.deg2rad(ang)

        # Push small values further out
        r = base_radius + (0.25 if val < small_thresh else 0.0)

        x = r * np.cos(ang_rad)
        y = r * np.sin(ang_rad)

        ax.text(x, y, f"{val:.1f}", ha="center", va="center", fontsize=7)

        # Leader line for small slices
        if val < small_thresh:
            x0 = 1.0 * np.cos(ang_rad)
            y0 = 1.0 * np.sin(ang_rad)
            ax.plot([x0, x], [y0, y], linewidth=0.5, color="black")


def plot_piecharts(
    df_car_pos,
    df_car_neg,
    palette,
    outpath="CAR_neighborhood_pies.png",
    annotate=True,  # <<< NEW
    fallback="#CCCCCC",
):
    fig_width_in = 6 / 2.54
    fig_height_in = 3 / 2.54

    colors_pos = get_colors_for_targets(
        df_car_pos["target"], palette, fallback=fallback
    )
    colors_neg = get_colors_for_targets(
        df_car_neg["target"], palette, fallback=fallback
    )

    wedge_kw = dict(edgecolor="white", linewidth=0.0)

    fig, axes = plt.subplots(
        1, 2, figsize=(fig_width_in, fig_height_in), constrained_layout=True
    )

    # --- CAR_pos ---
    wedges_pos, _ = axes[0].pie(
        df_car_pos["composition_percent"],
        labels=None,
        colors=colors_pos,
        startangle=90,
        wedgeprops=wedge_kw,
    )
    axes[0].set_title("CAR⁺ neighborhood", fontsize=7, pad=10)
    axes[0].axis("equal")

    if annotate:  # <<< conditional annotation
        annotate_pie(axes[0], wedges_pos, df_car_pos["neighbors_per_cell"])

    # --- CAR_neg ---
    wedges_neg, _ = axes[1].pie(
        df_car_neg["composition_percent"],
        labels=None,
        colors=colors_neg,
        startangle=90,
        wedgeprops=wedge_kw,
    )
    axes[1].set_title("CAR⁻ neighborhood", fontsize=7, pad=10)
    axes[1].axis("equal")

    if annotate:  # <<< conditional annotation
        annotate_pie(axes[1], wedges_neg, df_car_neg["neighbors_per_cell"])

    # --- shared legend ---
    fig.legend(
        wedges_pos,
        df_car_pos["target"],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
    )

    fig.savefig(outpath, dpi=600, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def plot_rl_interactions(
    adata,
    sample_id,
    comparison_dir,
    groupby="combined_celltype_and_CAR",
    source_labels=("CAR_pos", "CAR_neg"),
    target_labels=None,
    expr_prop=0.1,
    resource_name="consensus",
    uns_key="cpdb_diff_CAR_on_TME",
    interaction_diff_cutoff=0.1,
    pval_cutoff=0.05,
    width_cm=16.0,
    height_cm=10.0,
    dpi=600,
):

    # ----------------------------------------------------------
    # Filter interactions
    # ----------------------------------------------------------

    df = adata.uns["cpdb_res"].copy()

    df = df[df["target"].isin(target_labels)]

    df = df[df["source"].isin(source_labels)]

    interaction_cols = ["ligand", "receptor", "target"]

    pivot = df.pivot_table(
        index=interaction_cols,
        columns="source",
        values="lr_means",
        aggfunc="mean",
    )

    # ----------------------------------------------------------
    # Differential interactions
    # ----------------------------------------------------------

    gene_diff = pivot.isna().any(axis=1)

    value_diff = pivot.notna().all(axis=1) & (
        abs(pivot[source_labels[0]] - pivot[source_labels[1]]) > interaction_diff_cutoff
    )

    keep_interactions = pivot[gene_diff | value_diff].reset_index()

    df_filtered = df.merge(
        keep_interactions[interaction_cols],
        on=interaction_cols,
    )

    # ----------------------------------------------------------
    # Rename columns
    # ----------------------------------------------------------

    df_filtered = df_filtered.rename(
        columns={
            "lr_means": "expression magnitude",
            "cellphone_pvals": "−log10(p value)",
        }
    )

    # ----------------------------------------------------------
    # Store in AnnData
    # ----------------------------------------------------------

    adata.uns[uns_key] = df_filtered

    # ----------------------------------------------------------
    # Plot
    # ----------------------------------------------------------

    fig = li.pl.dotplot(
        adata=adata,
        colour="expression magnitude",
        size="−log10(p value)",
        inverse_size=True,
        source_labels=list(source_labels),
        target_labels=target_labels,
        filter_fun=lambda x: (x["−log10(p value)"] <= pval_cutoff),
        uns_key=uns_key,
        return_fig=True,
        size_range=(2, 6),
    )

    # ----------------------------------------------------------
    # Theme
    # ----------------------------------------------------------

    fig = fig + theme(
        figure_size=(
            width_cm / 2.54,
            height_cm / 2.54,
        ),
        text=element_text(size=7),
        axis_title_x=element_blank(),
        axis_title_y=element_blank(),
        axis_text_x=element_text(
            rotation=45,
            ha="right",
            va="top",
        ),
        axis_text_y=element_text(size=7),
        legend_title=element_text(size=6),
        legend_text=element_text(size=6),
        legend_position="right",
    )

    # ----------------------------------------------------------
    # Save
    # ----------------------------------------------------------

    figdir = os.path.join(
        comparison_dir,
        sample_id,
    )

    os.makedirs(figdir, exist_ok=True)

    outfile = os.path.join(
        figdir,
        f"{sample_id}_CPDB_CAR_dotplot.png",
    )

    fig.save(
        outfile,
        dpi=dpi,
        bbox_inches="tight",
    )

    print(f"Saved: {outfile}")

    return fig


def remove_small_clusters(adata, cluster_key, min_cells=150):
    """
    Remove clusters with fewer than `min_cells` cells from an AnnData object.
    """
    if cluster_key is None:
        raise ValueError("cluster_key must be provided.")

    print("Before filtering:")
    print("n_obs:", adata.n_obs)
    print(adata.obs[cluster_key].value_counts())

    # Count cells per cluster
    vc = adata.obs[cluster_key].value_counts()

    # Identify clusters to remove
    small_clusters = vc[vc < min_cells].index.astype(str).tolist()

    if small_clusters:
        labels_as_str = adata.obs[cluster_key].astype(str)
        adata = adata[~labels_as_str.isin(small_clusters)].copy()

        # Clean unused categories if categorical
        if pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
            adata.obs[cluster_key] = adata.obs[
                cluster_key
            ].cat.remove_unused_categories()

    print("\nRemoved clusters:", small_clusters if small_clusters else "None")

    print("\nAfter filtering:")
    print(adata.obs[cluster_key].value_counts())
    print("n_obs:", adata.n_obs)

    return adata


def add_celltype_and_CAR_column(adata, key):
    adata.obs[key] = np.where(
        adata.obs["Tcell_CAR_status"] != "non_T_cell",
        adata.obs["Tcell_CAR_status"],
        adata.obs["L2_fits_L1"],
    )
    return adata


def plot_neighbors_spatial(adata, cell_idx, title, crop_box):
    # initialize as strings
    adata.obs["_is_neighbor"] = "other"
    neighbors = adata.obsp["spatial_connectivities"][cell_idx].nonzero()[1]
    adata.obs.iloc[neighbors, adata.obs.columns.get_loc("_is_neighbor")] = "neighbor"

    adata.obs.iloc[cell_idx, adata.obs.columns.get_loc("_is_neighbor")] = "reference"

    # make categorical (optional but good practice)
    adata.obs["_is_neighbor"] = adata.obs["_is_neighbor"].astype("category")
    palette = {"reference": "#543937", "neighbor": "#508A88", "other": "#CEE5E7"}
    sc.pl.spatial(
        adata,
        color="_is_neighbor",
        title="",
        size=0.7,
        crop_coord=crop_box,
        palette=palette,
        legend_loc=None,
        frameon=False,
    )

    # cleanup
    adata.obs.drop(columns="_is_neighbor", inplace=True)


def get_top_interactions(adata, key, top=20):
    # Get z-score matrix
    z = adata.uns[f"{key}_nhood_enrichment"]["zscore"]
    cats = adata.obs[key].cat.categories

    # DataFrame
    z_df = pd.DataFrame(z, index=cats, columns=cats)

    # Long format
    z_long = z_df.stack().reset_index()
    z_long.columns = ["cluster1", "cluster2", "zscore"]

    # Only upper triangle, no self pairs
    z_long = z_long[z_long["cluster1"] < z_long["cluster2"]]

    # Keep only interactions where CAR_pos or CAR_neg is in either cluster
    mask = z_long["cluster1"].str.contains("CAR_pos|CAR_neg") | z_long[
        "cluster2"
    ].str.contains("CAR_pos|CAR_neg")
    z_filtered = z_long[mask]

    # Top N
    top_pairs = z_filtered.sort_values("zscore", ascending=False).head(top)
    return top_pairs


# --------------------------------------------------
# Plot neighborhood enrichment heatmap
# --------------------------------------------------
def plot_neighborhood_enrichment(
    adata,
    key,
    sample_id,
    output_dir=None,
    figsize=(10, 10),
    cmap="bwr",
    dpi=300,
    show=True,
):
    """
    Plot and optionally save Squidpy neighborhood
    enrichment heatmap.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    key : str
        Column used for neighborhood enrichment.

    sample_id : str
        Sample name for plot title.

    output_dir : str or Path, optional
        Directory to save figure.

    figsize : tuple
        Figure size.

    cmap : str
        Matplotlib colormap.

    dpi : int
        Figure DPI for saving.

    show : bool
        Whether to display the figure.
    """

    print(f"Plotting neighborhood enrichment for {sample_id}")

    enrichment = adata.uns[f"{key}_nhood_enrichment"]["zscore"]

    celltypes = adata.obs[key].cat.categories.tolist()

    # --------------------------------------------------
    # Normalize colormap around zero
    # --------------------------------------------------
    vmax = np.abs(enrichment).max()

    norm = mcolors.TwoSlopeNorm(
        vmin=-vmax,
        vcenter=0,
        vmax=vmax,
    )

    # --------------------------------------------------
    # Plot
    # --------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        enrichment,
        ax=ax,
        cmap=cmap,
        norm=norm,
        square=True,
        linewidths=0.5,
        linecolor="grey",
        xticklabels=celltypes,
        yticklabels=celltypes,
        cbar_kws={"label": "Z-score"},
    )
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(
        axis="both",
        which="both",
        length=0,
    )
    ax.set_title(f"{sample_id}\nNeighborhood Enrichment")

    ax.set_xlabel("Cell type")
    ax.set_ylabel("Cell type")

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)

    plt.tight_layout()

    # --------------------------------------------------
    # Save figure
    # --------------------------------------------------
    if output_dir is not None:

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        outfile = output_dir / f"{sample_id}/{key}_nhood_enrichment.png"

        fig.savefig(
            outfile,
            dpi=dpi,
            bbox_inches="tight",
        )

        print(f"Saved figure to:\n{outfile}")

    # --------------------------------------------------
    # Show / close
    # --------------------------------------------------
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_top_5_receptor_ligand_pairs(
    lrdata,
    sample_id,
    lr_pairs,
    crop_box,
    output_dir,
    dpi=600,
    show=True,
):
    """
    Plot top receptor-ligand spatial interaction pairs.

    Parameters
    ----------
    lrdata : AnnData
        Spatial ligand-receptor AnnData object.

    sample_id : str
        Sample identifier.

    lr_pairs : list
        List of ligand-receptor pairs to plot.

    crop_box : tuple
        Crop coordinates for spatial plot.

    output_dir : str or Path
        Base output directory.

    dpi : int
        Figure DPI.

    show : bool
        Whether to display figure.
    """

    n = len(lr_pairs)

    ncols = 3
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5 * ncols + 3, 5 * nrows),
    )

    axes = np.array(axes).flatten()

    # --------------------------------------------------
    # Plot LR pairs
    # --------------------------------------------------
    for i, pair in enumerate(lr_pairs):

        sc.pl.spatial(
            lrdata,
            layer="cats",
            color=[pair],
            size=0.8,
            crop_coord=crop_box,
            frameon=False,
            cmap="coolwarm",
            ax=axes[i],
            show=False,
            colorbar_loc=None,
        )

        axes[i].set_title(
            pair,
            fontsize=14,
        )

    # --------------------------------------------------
    # Hide unused axes
    # --------------------------------------------------
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    # --------------------------------------------------
    # Custom legend
    # --------------------------------------------------
    legend_elements = [
        Patch(facecolor="#b40426", label="High–High"),
        Patch(facecolor="#e0e0e0", label="No association"),
        Patch(facecolor="#3b4cc0", label="High–Low/Low-High"),
    ]

    fig.legend(
        handles=legend_elements,
        loc="center right",
        frameon=False,
        title="Local LR category",
        fontsize=14,
        title_fontsize=14,
    )

    # --------------------------------------------------
    # Layout
    # --------------------------------------------------
    plt.subplots_adjust(
        wspace=0.1,
        hspace=0.02,
        right=0.85,
    )

    # --------------------------------------------------
    # Save figure
    # --------------------------------------------------
    output_dir = Path(output_dir) / sample_id / "local_niches"

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    outfile = output_dir / f"{sample_id}_top_lr_pairs.png"

    fig.savefig(
        outfile,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
    )

    print(f"Saved figure to:\n{outfile}")

    # --------------------------------------------------
    # Show / close
    # --------------------------------------------------
    if show:
        plt.show()
    else:
        plt.close(fig)
