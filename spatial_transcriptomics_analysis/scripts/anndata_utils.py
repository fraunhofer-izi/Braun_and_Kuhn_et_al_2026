from __future__ import annotations
import scanpy as sc
import HERGAST
import numpy as np
from typing import Iterable, Optional, Sequence, Union, Dict, List
import pandas as pd
from pathlib import Path
import os
import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
logger = logging.getLogger(__name__)


def attach_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    adata.var_names_make_unique()
    adata = adata[:, ~adata.var_names.str.startswith("DEPRECATED")].copy()
    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito", "ribo"], inplace=True)
    return adata

def annotate_with_panglao(
    adata: sc.AnnData,
    panglao_file: Union[str, Path] = "/mnt/ribolution/user_worktmp/christina.kuhn/work/references/PanglaoDB_markers_27_Mar_2020.tsv",
    species: str = "Hs",
    organs_of_interest: Optional[Sequence[str]] = None,
    specific_celltypes: Optional[Sequence[str]] = None,
    score_threshold: float = 0.1,
    score_prefix: str = "panglao",
    inplace: bool = True,
) -> sc.AnnData:
    """
    Annotate cells with PanglaoDB marker scores and choose the best-matching cell type.

    Parameters
    ----------
    adata
        AnnData object to annotate.
    panglao_file
        Path to PanglaoDB TSV (e.g. PanglaoDB_markers_27_Mar_2020.tsv).
    species
        Species string used in PanglaoDB (e.g. "Hs", "Mm"). Substring match.
    organs_of_interest
        If provided, restrict markers to these organs.
    specific_celltypes
        If provided, restrict markers to these cell types.
    score_threshold
        Minimum max-score required to assign a cell type; otherwise "Undetermined".
    score_prefix
        Prefix for created columns in `adata.obs` (e.g. `<prefix>_score_*`).
    inplace
        If False, operate on a copy and return it.

    Returns
    -------
    AnnData
        The (modified) AnnData with:
        - obs[f"{score_prefix}_score_<celltype>"]  : per-cell marker scores
        - obs[f"{score_prefix}_celltype"]          : top-scoring label (raw)
        - obs[f"{score_prefix}_scoremax"]          : max score value
        - obs[f"{score_prefix}_celltype_filtered"] : label after thresholding
    """
    if not inplace:
        adata = adata.copy()

    panglao_path = Path(panglao_file)
    if not panglao_path.exists():
        raise FileNotFoundError(f"PanglaoDB file not found: {panglao_path}")

    # Load & validate input table
    df = pd.read_csv(panglao_path, sep="\t")
    required_cols = {"species", "organ", "cell type", "official gene symbol"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(
            f"PanglaoDB table missing required columns: {', '.join(sorted(missing))}"
        )

    # Filter by species (substring match to allow composite values)
    df = df[df["species"].astype(str).str.contains(species, na=False, regex=False)].copy()

    # Optional filters
    if organs_of_interest is not None:
        organs_set = set(organs_of_interest)
        df = df[df["organ"].isin(organs_set)].copy()

    # Build mapping: cell type -> unique, sorted gene symbols
    df = df.dropna(subset=["official gene symbol", "cell type"]).copy()
    df["official gene symbol"] = df["official gene symbol"].astype(str)
    df["cell type"] = df["cell type"].astype(str)

    celltype_to_genes: Dict[str, List[str]] = (
        df.groupby("cell type")["official gene symbol"]
        .apply(lambda s: sorted(set(s.tolist())))
        .to_dict()
    )

    if specific_celltypes is not None:
        keep = set(specific_celltypes)
        celltype_to_genes = {ct: genes for ct, genes in celltype_to_genes.items() if ct in keep}

    # Intersect with genes present in the AnnData
    valid_genes = set(map(str, adata.var_names))
    celltype_to_genes = {
        ct: [g for g in genes if g in valid_genes]
        for ct, genes in celltype_to_genes.items()
    }
    # Drop any empty sets to avoid score_genes errors
    celltype_to_genes = {ct: genes for ct, genes in celltype_to_genes.items() if len(genes) > 0}

    if not celltype_to_genes:
        raise ValueError(
            "No Panglao markers overlap with `adata.var_names` after filtering. "
            "Check `species`, `organs_of_interest`, `specific_celltypes`, or gene IDs."
        )

    # Helper to make safe score names
    def _slugify(s: str) -> str:
        return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_")

    # Score genes for each cell type
    score_name_by_ct: Dict[str, str] = {}
    for ct, genes in celltype_to_genes.items():
        score_name = f"{score_prefix}_score_{_slugify(ct)}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_name)
        score_name_by_ct[ct] = score_name

    # Collect the created score columns (guard in case some already existed)
    score_keys = [nm for nm in score_name_by_ct.values() if nm in adata.obs.columns]
    if not score_keys:
        raise RuntimeError("No score columns were created in `adata.obs`.")

    scores = adata.obs[score_keys]
    adata.obsm["celltype_scores"] = scores.copy()

    # Compute top label and max value with pandas (keeps index aligned)
    top_score_col = scores.idxmax(axis=1)
    max_score_val = scores.max(axis=1)

    score_to_celltype = {v: k for k, v in score_name_by_ct.items()}
    top_score_col_clean = top_score_col.map(score_to_celltype)

    # Apply cutoff (e.g., 0.1) and store results
    th = float(score_threshold)  # e.g., 0.1
    labels = (
        top_score_col_clean.astype(object)
        .where(max_score_val >= th, other="Undetermined")
        .astype("category")
    )

    adata.obs[f"{score_prefix}_celltype"] = labels
    adata.obs[f"{score_prefix}_scoremax"] = max_score_val.astype(float)

    return adata

def annotate_with_panglao_aucell(
    adata: sc.AnnData,
    panglao_file: Union[str, Path] = "/mnt/ribolution/user_worktmp/christina.kuhn/work/references/PanglaoDB_markers_27_Mar_2020.tsv",
    species: str = "Hs",
    organs_of_interest: Optional[Sequence[str]] = None,
    specific_celltypes: Optional[Sequence[str]] = None,
    score_threshold: float = 0.1,            # AUCell is in [0,1]; 0.05–0.1 are common starting points
    score_prefix: str = "panglao",
    min_markers_per_set: int = 5,            # drop very small sets (AUCell is fine, but tiny sets are noisy)
    use_raw: bool = False,                    # pass through to decoupler (rank-based, raw not required)
    inplace: bool = True,
) -> sc.AnnData:
    """
    Annotate cells with PanglaoDB using AUCell (rank-based, robust to dropout).

    Produces:
      - obs[f"{score_prefix}_score_<celltype>"]  : per-cell AUCell scores (0..1)
      - obs[f"{score_prefix}_celltype"]          : top-scoring label (raw)
      - obs[f"{score_prefix}_scoremax"]          : max AUCell score value
      - obs[f"{score_prefix}_celltype_filtered"] : label after thresholding
    """
    if not inplace:
        adata = adata.copy()

    panglao_path = Path(panglao_file)
    if not panglao_path.exists():
        raise FileNotFoundError(f"PanglaoDB file not found: {panglao_path}")

    # --- Load & filter PanglaoDB
    df = pd.read_csv(panglao_path, sep="\t")
    required_cols = {"species", "organ", "cell type", "official gene symbol"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"PanglaoDB table missing required columns: {', '.join(sorted(missing))}")

    # species filter (substring to allow composite notations)
    df = df[df["species"].astype(str).str.contains(species, na=False, regex=False)].copy()

    if organs_of_interest is not None:
        organs_set = set(organs_of_interest)
        df = df[df["organ"].isin(organs_set)].copy()

    df = df.dropna(subset=["official gene symbol", "cell type"]).copy()
    df["official gene symbol"] = df["official gene symbol"].astype(str)
    df["cell type"] = df["cell type"].astype(str)

    # Build mapping: cell type -> unique, sorted gene symbols
    celltype_to_genes: Dict[str, List[str]] = (
        df.groupby("cell type")["official gene symbol"]
        .apply(lambda s: sorted(set(s.tolist())))
        .to_dict()
    )

    if specific_celltypes is not None:
        keep = set(specific_celltypes)
        celltype_to_genes = {ct: genes for ct, genes in celltype_to_genes.items() if ct in keep}

    # Intersect with genes present in AnnData
    valid_genes = set(map(str, adata.var_names))
    celltype_to_genes = {
        ct: [g for g in genes if g in valid_genes]
        for ct, genes in celltype_to_genes.items()
    }

    # Drop empty or tiny sets
    celltype_to_genes = {
        ct: genes for ct, genes in celltype_to_genes.items() if len(genes) >= min_markers_per_set
    }

    if not celltype_to_genes:
        raise ValueError(
            "No Panglao markers overlap (or pass min_markers_per_set) with `adata.var_names` after filtering."
        )

    # --- AUCell with decoupler expects a long-format network: columns ['source','target','weight']
    net = pd.DataFrame(
        [(ct, g, 1.0) for ct, glist in celltype_to_genes.items() for g in glist],
        columns=["source", "target", "weight"],
    )

    # Run AUCell (rank-based enrichment per cell)
    # Results: adata.obsm['aucell_estimate'] with columns per 'source' (cell type)
    dc.mt.aucell(
        adata,
        net,
        tmin=5
    )
    aucell = adata.obsm["score_aucell"].copy()

    # --- Write per-celltype scores into obs with your naming scheme
    def _slugify(s: str) -> str:
        return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_")

    score_cols = {}
    for ct in aucell.columns:
        col = f"{score_prefix}_score_{_slugify(ct)}"
        adata.obs[col] = aucell[ct].astype(float).values
        score_cols[ct] = col

    # Collect back the created columns in the same order as aucell.columns
    obs_scores = adata.obs[[score_cols[ct] for ct in aucell.columns]]

    # Top label and max score
    top_score_col = obs_scores.idxmax(axis=1)
    max_score_val = obs_scores.max(axis=1)

    # map obs column name -> original cell type
    inv_map = {v: k for k, v in score_cols.items()}
    top_ct = top_score_col.map(inv_map)

    adata.obs[f"{score_prefix}_celltype"] = top_ct.astype("category")
    adata.obs[f"{score_prefix}_scoremax"] = max_score_val.astype(float)

    # Threshold and filtered label
    th = float(score_threshold)
    filt_label = top_ct.where(max_score_val >= th, other="Undetermined").astype("category")
    adata.obs[f"{score_prefix}_celltype_filtered"] = filt_label

    return adata


import numpy as np
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
import numpy as np



def find_first_valley(scores,
                      min_points: int = 50,
                      grid_size: int = 800):
    """
    Findet das erste sinnvolle Tal (lokales Minimum) zwischen den beiden
    links-liegenden Hauptpeaks der Dichteverteilung.
    Funktioniert auch, wenn der zweite Peak höher ist als der erste.

    Parameters
    ----------
    scores : array-like
        1D array mit AUCell-Scores.
    min_points : int
        Minimum number of data points needed for KDE.
    grid_size : int
        Anzahl der Stützstellen für die KDE.
    """

    scores = np.asarray(scores)
    scores = scores[np.isfinite(scores)]
    if scores.size < min_points:
        return None

    # KDE auf vollständiger Verteilung
    kde = gaussian_kde(scores)
    xs = np.linspace(scores.min(), scores.max(), grid_size)
    ys = kde(xs)

    # Lokale Maxima (Peaks) und Minima (Täler)
    maxima_idx = argrelextrema(ys, np.greater)[0]
    minima_idx = argrelextrema(ys, np.less)[0]

    if maxima_idx.size == 0 or minima_idx.size == 0:
        return None

    # globalen Max-Peak sicherheitshalber hinzufügen, falls argrelextrema ihn auslässt
    global_max_idx = int(np.argmax(ys))
    if global_max_idx not in maxima_idx:
        maxima_idx = np.sort(np.append(maxima_idx, global_max_idx))

    # Peaks nach x sortieren und kleine Peaks rausfiltern
    maxima_idx = np.array(sorted(maxima_idx))

    # brauchen mindestens 2 "echte" Peaks
    if maxima_idx.size < 2:
        return None

    # nimm die ersten beiden Peaks (links-liegend)
    first_peak_idx = maxima_idx[0]
    second_peak_idx = maxima_idx[1]

    # alle Minima zwischen erstem und zweitem Peak
    between = minima_idx[(minima_idx > first_peak_idx) & (minima_idx < second_peak_idx)]
    if between.size == 0:
        return None

    # tiefstes Minimum zwischen den beiden Peaks
    valley_idx = between[np.argmin(ys[between])]
    valley_x = float(xs[valley_idx])
    return valley_x



def annotate_with_custom_aucell(
    adata: sc.AnnData,
    json_file: Union[str, Path],
    score_threshold: float = 0.1,
    score_prefix: str = "custom",
    min_markers_per_set: int = 5,
    use_raw: bool = False,
    inplace: bool = True,
) -> sc.AnnData:

    import json
    import seaborn as sns
    import matplotlib.pyplot as plt

    if not inplace:
        adata = adata.copy()

    json_path = Path(json_file)
    if not json_path.exists():
        raise FileNotFoundError(f"Marker JSON file not found: {json_path}")

    # -------------------------
    # Load markers
    # -------------------------
    with open(json_path, "r", encoding="utf-8") as f:
        celltype_to_genes_raw = json.load(f)

    valid_genes = set(map(str, adata.var_names))

    # Track missing genes
    missing_gene_info = []

    celltype_to_genes = {}
    for ct, genes in celltype_to_genes_raw.items():
        genes = list(map(str, genes))
        present = [g for g in genes if g in valid_genes]
        missing = [g for g in genes if g not in valid_genes]

        if missing:
            missing_gene_info.append({"celltype": ct, "missing_genes": ",".join(missing)})

        if len(present) >= min_markers_per_set:
            celltype_to_genes[ct] = present

    if not celltype_to_genes:
        raise ValueError("No overlapping marker genes.")

    # Store missing gene information
    adata.uns[f"{score_prefix}_missing_genes"] = pd.DataFrame(missing_gene_info)

    # -------------------------
    # Create network table
    # -------------------------
    net = pd.DataFrame(
        [(ct, g, 1.0) for ct, genes in celltype_to_genes.items() for g in genes],
        columns=["source", "target", "weight"],
    )

    # -------------------------
    # Run AUCell
    # -------------------------
    dc.mt.aucell(adata, net, tmin=5)
    aucell = adata.obsm["score_aucell"].copy()

    # -------------------------
    # Write AUCell scores into obs
    # -------------------------
    def _slugify(s: str) -> str:
        return "".join(ch if ch.isalnum() else "_" for ch in s).strip("_")

    score_cols = {}
    for ct in aucell.columns:
        col = f"{score_prefix}_score_{_slugify(ct)}"
        adata.obs[col] = aucell[ct].astype(float).values
        score_cols[ct] = col

    obs_scores = adata.obs[[score_cols[ct] for ct in aucell.columns]]

    # -------------------------
    # Compute per-celltype quantile thresholds
    # -------------------------
    thresholds = {}

    # For plotting
    plot_dict = {}

    for ct in aucell.columns:
        col = score_cols[ct]
        scores = obs_scores[col].values

        # Versuche automatischen Talpunkt
        valley = find_first_valley(scores)
        if valley is None:
            valley = 0.1
        thresholds[ct] = valley

        # -------------------------
        # Plot the distribution
        # -------------------------
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.histplot(scores, ax=ax, kde=True, bins=50)

        # draw cutoff only if valley was found
        if valley is not None:
            ax.axvline(valley, color="red", linestyle="--", linewidth=2)
            title_line = f"valley cutoff = {valley:.4f}"
        else:
            title_line = "valley cutoff = None"

        ax.set_title(f"AUCell Score Distribution – {ct}\n{title_line}")
        ax.set_xlabel("AUCell score")
        ax.set_ylabel("Count")
        plot_dict[ct] = fig

    # Save plots
    adata.uns[f"{score_prefix}_aucell_plots"] = plot_dict

    # Show plots inline
    for fig in plot_dict.values():
        fig.show()

    print("Computed quantile thresholds:\n" +
          "\n".join([f"  {ct}: {thr:.4f}" for ct, thr in thresholds.items()]))

    # -------------------------
    # Assign celltypes based on thresholds
    # -------------------------
    celltype_filtered = []
    celltypes_passing_all = []

    for idx in obs_scores.index:
        row = obs_scores.loc[idx]

        # --- which CTs pass threshold?
        passing = []
        for ct in aucell.columns:
            score = row[score_cols[ct]]
            if score >= thresholds[ct]:
                passing.append(ct)

        # store all passing (new requirement)
        celltypes_passing_all.append(",".join(passing) if passing else "")

        if len(passing) == 0:
            celltype_filtered.append("Undetermined")
            continue

        # choose highest scoring among passing
        best = max(passing, key=lambda ct: row[score_cols[ct]])
        celltype_filtered.append(best)

    adata.obs[f"{score_prefix}_celltype_filtered"] = pd.Categorical(celltype_filtered)
    adata.obs[f"{score_prefix}_celltypes_passing"] = celltypes_passing_all

    return adata


def annotate_with_custom_aucell_threshold(
    adata: sc.AnnData,
    json_file: Union[str, Path],
    score_threshold: float = 0.1,
    score_prefix: str = "custom",
    min_markers_per_set: int = 5,
    inplace: bool = True,
) -> sc.AnnData:

    import json
    import seaborn as sns
    import matplotlib.pyplot as plt

    if not inplace:
        adata = adata.copy()

    json_path = Path(json_file)
    if not json_path.exists():
        raise FileNotFoundError(f"Marker JSON file not found: {json_path}")

    # -------------------------
    # Load markers
    # -------------------------
    with open(json_path, "r", encoding="utf-8") as f:
        celltype_to_genes_raw = json.load(f)

    valid_genes = set(map(str, adata.var_names))

    missing_gene_info = []
    celltype_to_genes = {}

    for ct, genes in celltype_to_genes_raw.items():
        genes = list(map(str, genes))
        present = [g for g in genes if g in valid_genes]
        missing = [g for g in genes if g not in valid_genes]

        if missing:
            missing_gene_info.append({"celltype": ct, "missing_genes": ",".join(missing)})

        if len(present) >= min_markers_per_set:
            celltype_to_genes[ct] = present

    if not celltype_to_genes:
        raise ValueError("No marker sets contain ≥ min_markers_per_set present genes.")

    adata.uns[f"{score_prefix}_missing_genes"] = pd.DataFrame(missing_gene_info)

    # -------------------------
    # Build AUCell network
    # -------------------------
    net = pd.DataFrame(
        [(ct, g, 1.0) for ct, genes in celltype_to_genes.items() for g in genes],
        columns=["source", "target", "weight"],
    )

    # -------------------------
    # Run AUCell
    # -------------------------
    dc.mt.aucell(adata, net, tmin=5)

    # rename Aucell score
    aucell = adata.obsm["score_aucell"].copy()
    adata.obsm[f"{score_prefix}_aucell"] = adata.obsm["score_aucell"]

    # -------------------------
    # Compute "filtered" and "passing" annotations
    # -------------------------
    celltype_filtered = []
    celltypes_passing_all = []

    for idx in aucell.index:
        row = aucell.loc[idx]

        # CTs above threshold
        passing = [ct for ct in aucell.columns if row[ct] >= score_threshold]

        celltypes_passing_all.append(",".join(passing) if passing else "")

        if len(passing) == 0:
            celltype_filtered.append("Undetermined")
        else:
            best = max(passing, key=lambda ct: row[ct])
            celltype_filtered.append(best)

    # Store final annotations (no score columns!)
    adata.obs[f"{score_prefix}_celltype_filtered"] = pd.Categorical(celltype_filtered)
    adata.obs[f"{score_prefix}_celltypes_passing"] = celltypes_passing_all

    return adata

def check_l2_matches_l1(
    adata,
    json_file: str,
    l2_col: str,
    l1_passing_col: str,
    output_col: str = "max_celltype_fits_l1",
    l1_from_l2_col: str = "l1_from_l2",
    fits_boolean_col: str = "l2_l1_consistent",
):
    """
    Check whether the L2 max celltype maps (via L1→L2 JSON) to an L1 celltype
    that is included in the L1 passing list.
    """

    import json
    import numpy as np
    import pandas as pd

    # ----------------------
    # Load L1 → [L2 list] mapping
    # ----------------------
    with open(json_file, "r") as f:
        l1_to_l2 = json.load(f)

    # ----------------------
    # Reverse mapping: L2 → L1
    # ----------------------
    l2_to_l1 = {}
    for l1, l2_list in l1_to_l2.items():
        for l2 in l2_list:
            l2_to_l1[str(l2)] = str(l1)

    # ----------------------
    # Robustly parse the L1 passing column
    # Handles: list, empty list, "", NaN, string "A,B,C"
    # ----------------------
    def parse_list(x):
        if isinstance(x, list):
            # already a list → clean it
            return [str(i).strip() for i in x]
        if pd.isna(x) or x is None:
            return []
        if isinstance(x, str):
            if x.strip() == "":
                return []
            # split comma-separated string
            return [i.strip() for i in x.split(",")]
        # unexpected type → fallback to empty
        return []

    # apply without pandas trying to hash lists
    l1_passing_lists = adata.obs[l1_passing_col].map(parse_list)

    mapped_l1_results = []
    fits_boolean_results = []
    output_results = []

    # ----------------------
    # Iterate through cells
    # ----------------------
    for idx, l2_value in adata.obs[l2_col].items():

        l2_value_str = str(l2_value)

        # Map L2 → L1 (None if missing)
        mapped_l1 = l2_to_l1.get(l2_value_str, None)
        mapped_l1_results.append(mapped_l1)

        passing_l1 = l1_passing_lists.loc[idx]

        # Boolean: does mapped L1 appear in passing list?
        fits = (mapped_l1 is not None) and (mapped_l1 in passing_l1)
        fits_boolean_results.append(fits)

        # Output L2 if fits, else NaN
        output_results.append(l2_value if fits else np.nan)

    # ----------------------
    # Save results to AnnData
    # ----------------------
    adata.obs[l1_from_l2_col] = mapped_l1_results
    adata.obs[fits_boolean_col] = fits_boolean_results
    adata.obs[output_col] = output_results

    return adata

def cluster_based_on_gex(adata: sc.AnnData) -> sc.AnnData:
    sc.settings.n_jobs = -1
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata_scaled = adata.copy()
    sc.pp.highly_variable_genes(adata_scaled, n_top_genes=1500, flavor="seurat_v3", inplace=True)
    adata_scaled = adata_scaled[:, adata_scaled.var.highly_variable]
    sc.pp.scale(adata_scaled)

    for method in ["pca", "neighbors", "umap", "leiden"]:
        if method == "pca":
            sc.tl.pca(adata_scaled, n_comps=100)  # compute more PCs than needed
            vr = adata_scaled.uns['pca']['variance_ratio']
            dvr = np.diff(vr)  # drop in variance explained
            # first PC index where the drop gets small and stays small
            n_pcs = int(np.argmax(dvr < 0.001) + 1) or 30
            n_pcs = max(20, min(n_pcs, 60))  # enforce bounds
            adata_scaled.uns['n_pcs_opt'] = n_pcs  # store for later
        elif method == "neighbors":
            sc.pp.neighbors(adata_scaled, n_neighbors=15,
                            n_pcs=int(adata_scaled.uns.get('n_pcs_opt', 50)))
        elif method == "umap":
            sc.tl.umap(adata_scaled)
        elif method == "leiden":
            resolutions = np.arange(0.1, 1.1, 0.1)
            for res in resolutions:
                key = f"leiden_{res:.1f}".replace(".", "_")  # e.g., leiden_0_1
                sc.tl.leiden(adata_scaled, resolution=res, key_added=key, random_state=2024)
    # Copy embeddings & clustering result back
    adata.obsm["X_pca_GEX"] = adata_scaled.obsm["X_pca"]
    adata.obsm["X_umap_GEX"] = adata_scaled.obsm["X_umap"]

    adata.uns["pca_GEX"] = adata_scaled.uns["pca"]
    adata.uns["neighbors_GEX"] = adata_scaled.uns["neighbors"]
    adata.uns["umap_GEX"] = adata_scaled.uns["umap"]

    adata.obsp["connectivities_GEX"] = adata_scaled.obsp["connectivities"]
    adata.obsp["distances_GEX"] = adata_scaled.obsp["distances"]

    resolutions = np.arange(0.1, 1.1, 0.1)
    for res in resolutions:
        key = f"leiden_{res:.1f}".replace(".", "_")
        adata.obs[key] = adata_scaled.obs[key]
    adata.uns['n_pcs_opt'] = adata_scaled.uns['n_pcs_opt']
    return adata

def clustering_with_hergast(adata: sc.AnnData) -> sc.AnnData:
    tmp = adata.copy()
    if tmp.raw is not None:
        print("Detected .raw — using raw counts from .raw.X for HERGAST preprocessing.")
        tmp.X = tmp.raw.X.copy()

    sc.pp.normalize_total(tmp, target_sum=1, exclude_highly_expressed=True)
    sc.pp.scale(tmp)
    sc.pp.pca(tmp, n_comps=100)
    HERGAST.utils.Cal_Spatial_Net(tmp)
    HERGAST.utils.Cal_Expression_Net(tmp, dim_reduce='PCA')
    train_HERGAST = HERGAST.Train_HERGAST(tmp, batch_data=True, num_batch_x_y=(7,7), spatial_net_arg={'verbose':False},
                                      exp_net_arg={'verbose':False},dim_reduction='PCA',device_idx=0)
    train_HERGAST.train_HERGAST(n_epochs=200)
    sc.pp.neighbors(tmp, use_rep='HERGAST')
    sc.tl.umap(tmp)
    sc.tl.leiden(tmp, random_state=2024, resolution=0.3, key_added='leiden_HERGAST_0_3')
    sc.tl.leiden(tmp, random_state=2024, resolution=0.4, key_added='leiden_HERGAST_0_4')
    sc.tl.leiden(tmp, random_state=2024, resolution=0.5, key_added='leiden_HERGAST_0_5')

    adata.obsp["connectivities_HERGAST"] = tmp.obsp["connectivities"].copy()
    if "distances" in tmp.obsp:
        adata.obsp["distances_HERGAST"] = tmp.obsp["distances"].copy()
    adata.uns["neighbors_HERGAST"] = tmp.uns.get("neighbors", {}).copy()
    adata.obs["leiden_HERGAST_0_3"] = tmp.obs["leiden_HERGAST_0_3"].copy()
    adata.obs["leiden_HERGAST_0_4"] = tmp.obs["leiden_HERGAST_0_4"].copy()
    adata.obs["leiden_HERGAST_0_5"] = tmp.obs["leiden_HERGAST_0_5"].copy()
    adata.obsm["X_umap_HERGAST"] = tmp.obsm["X_umap"].copy()
    adata.obsm["X_pca_HERGAST"] = tmp.obsm["X_pca"].copy()
    adata.uns['HERGAST'] = {
        'n_epochs': 200,
        'batch_data': True,
        'num_batch_x_y':  [7, 7],
        'spatial_net_arg': {'verbose': False},
        'exp_net_arg': {'verbose': False},
        'dim_reduction': 'PCA',
        'device_idx': 0
    }
    return adata

def extract_metric_df(adatas: dict, metric: str) -> pd.DataFrame:
    """
    Combines a given metric across multiple AnnData objects into one DataFrame.
    """
    data = []
    for sample_name, adata in adatas.items():
        if metric not in adata.obs.columns:
            raise ValueError(f"Metric '{metric}' not found in '{sample_name}'")
        df = pd.DataFrame({
            "value": adata.obs[metric],
            "sample": sample_name
        })
        data.append(df)
    return pd.concat(data, ignore_index=True)


def plot_metric_boxplot(df: pd.DataFrame, metric: str, output_dir: str):
    """
    Creates and saves a boxplot for a given metric DataFrame.
    """
    plt.figure(figsize=(8, 5))
    sns.boxplot(data=df, x="sample", y="value")
    plt.ylabel(metric.replace("_", " ").capitalize())
    plt.title(f"{metric.replace('_', ' ').capitalize()} Across Samples")
    sns.despine(top=True, right=True)
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{metric}_boxplot_comparison.png")
    plt.savefig(out_path, dpi=300)
    plt.show()
    plt.close()
    print(f"✅ Saved: {out_path}")