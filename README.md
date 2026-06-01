# Braun and Kuhn et al. 2026 – Spatial Transcriptomics Analysis

This repository contains the code used to reproduce the spatial transcriptomics analysis presented in:

> Braun and Kuhn et al. – *Multimodal Profiling Reveals Actionable Vulnerabilities in CAR T–Derived T-Cell Lymphoma*
> (manuscript not yet released)

---

# Repository Structure

```text
.
├── assets/                     # Marker gene sets, color maps, mappings
├── config/                     # Sample-specific YAML configuration files
├── data/                       # Input data directories
├── notebooks/                  # **Figure generation notebook**
├── results/                    # Intermediate and final outputs
├── scripts/                    # Main analysis pipeline scripts
├── spatial_transcriptomics_analysis/
│   └── scripts/                # Reusable helper modules
├── environment.yml             # Conda environment
├── pyproject.toml              # Python package configuration
├── Makefile                    # Pipeline execution commands
└── README.md
```

---

# Reproducibility

## 1. Clone Repository

```bash
git clone <repository-url>
cd Braun_and_Kuhn_et_al_2026
```

---

## 2. Create Environment and install local package

```bash
make conda_env_create
```
This command will:

Create or update the Mamba environment
Install the local `spatial_transcriptomics_analysis package in editable mode

### Requirements

- conda/miniforge installed 
- mamba installed
- GNU make installed

tested with:
mamba 1.4.1
conda 23.1.0
GNU make 4.3
---

## 3. Download data

Download the GEO dataset:

```bash
mkdir -p data/spaceranger/
cd data/spaceranger/

wget -O GSE317410_RAW.tar \
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE317410&format=file"

tar -xvf GSE317410_RAW.tar
```

After extraction, organize the processed SpaceRanger output directory as:

```text
data/spaceranger/Duodenum_run3/outs
```

The pipeline assumes this directory structure for downstream processing.


## 4. Run Pipeline 

The full analysis pipeline can be executed using:

```bash
make all
```

The pipeline consists of the following stages:

| Step                    | Script                 | Description                       |
| ----------------------- | ---------------------- | --------------------------------- |
| Read SpaceRanger output | `01_read_spaceranger/` | Create SpatialData objects        |
| Quality control         | `02_qc/`               | Cell and gene filtering           |
| Clustering              | `03_cluster/`          | Spatial clustering and embedding  |
| Annotation              | `04_annotation/`       | Cell type annotation using AUCell |

Intermediate results are stored in:

```text
results/
├── zarr_before_qc/
├── zarr_after_qc/
├── zarr_after_cluster/
├── zarr_after_anno/
└── figures/
```

---

# 4. Publication Figures

Publication figures are generated using:

```text
notebooks/PlotsForPublication.ipynb
```

The notebook reproduces:

* neighborhood enrichment analyses
* ligand–receptor interaction plots
* spatial niche visualizations
* marker gene plots
* publication-ready summary figures

Generated figures are stored in:

```text
results/publication_plots/
```

---

# Configuration Files

Sample-specific settings are stored in:

```text
config/
├── Duodenum_run3.yaml
└── Skin.yaml
```

These files define:

* sample names
* clustering parameters
* annotation settings

---

# Assets

The `assets/` directory contains curated resources used for annotation and plotting:

* marker gene sets
* cell type mapping dictionaries
* color palettes

---

# Software Versions

The analysis was developed and tested using:

* Python 3.11
* Scanpy
* Squidpy
* SpatialData
* LIANA

Exact package versions of all packages are defined in:

```text
environment.yml
```

---

# Notes

* The repository uses editable installation via `pyproject.toml`.

---

# Citation

If you use this repository, please cite:

> Braun and Kuhn et al. – *Multimodal Profiling Reveals Actionable Vulnerabilities in CAR T–Derived T-Cell Lymphoma*

---

# Contact

For questions regarding the analysis pipeline, please open an issue in the repository.

---

# License

Copyright 2025 Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.

Licensed under the GPL-3.0. You may obtain a copy of the License in the LICENSE file.
