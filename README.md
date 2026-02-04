# Braun and Kuhn et al. 2026 – Spatial Transcriptomics Analysis

This repository contains the code used to produce the results for the spatial transcriptomics analysis:

> Braun and Kuhn et al. – Multimodal Profiling Reveals Actionable Vulnerabilities in CAR T–Derived T-Cell Lymphoma 
> (manuscript not yet released)

---

## Overview

This project provides the analysis pipeline for **spatial transcriptomics data** generated using the **10x Visium HD** platform.  
The workflows cover preprocessing of SpaceRanger output, clustering of cells, annotation of celltypes, and downstream analysis required to reproduce the study results.

---

## Analysis Structure

### Visium HD Data Analysis

The analysis of Visium HD samples was conducted using **Python**.
The pipeline consists of several modular scripts (Makefile) and two notebooks for downstream analysis and figure creation (Jupyter Notebooks):

| Step | Script | Description |
|------|--------|-------------|
| Read SpaceRanger output | `read_spaceRangerv4.py` | Reads SpaceRanger v4 output including segmented cells into a SpatialData object |
| Quality Control | `spatial_data_qc.py` | Does QC control steps of cells/genes |
| Clustering | `norm_and_cluster.py` | Performs spatially-aware clustering of cells |
| Annoation | `annotation_with_aucell.py` | Annotates cell types using AuCell |

| Notebook | Script | Description |
|------|--------|-------------|
| Spatial Regions and Local Niches | `01_SpatialRegions_and_LocalNiches.ipyb` | Plots generated from Spatial Regions and Local Niches Analysis  |
| CAR T microenvironment | `02_CART_microenvironemnt.ipyb` | Plots generated from CAR T microenvironment analysis (Neighboorhood, DEGs, Receptor-Ligand) |

---

## Reproducibility & Environment Setup

We use **Conda** for dependency management.

### Requirements

- Install **Miniconda** or **Anaconda**

### Create the Environment

A predefined Conda environment file is included in the repository.  
You can recreate the environment by running:

```bash
make conda_env_create ENV_NAME=spatialenv
```

### Run 

```
make all ENV_NAME=spatialenv
```

## Coding Guidelines

No strict style guide was enforced during development.  
However, code formatting consistency is maintained using:

- **Black** (Python code formatter)

---