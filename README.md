# Braun et al. 2026 – Spatial Transcriptomics Analysis

This repository contains the code used to produce the results for the spatial transcriptomics analysis:

> Braun et al. – Multimodal Profiling Reveals Actionable Vulnerabilities in CAR T–Derived T-Cell Lymphoma 
> (manuscript not yet released)

---

## Overview

This project provides analysis pipelines for **spatial transcriptomics data** generated using the **10x Genomics Visium** and **Visium HD** platforms.  
The workflows cover preprocessing, annotation, and downstream analysis required to reproduce the study results.

---

## Analysis Structure

### Visium Data Analysis (Lower Resolution)

The analysis of standard Visium samples was conducted in **R**.

Please refer to the dedicated README inside the corresponding Visium analysis directory for:

- Installation instructions  
- Required dependencies  
- Usage examples  
- Output descriptions  

---

### Visium HD Data Analysis (Higher Resolution)

The analysis of Visium HD samples was conducted using a combination of **Python** and **R**.  
The pipeline consists of several modular scripts:

| Step | Script | Description |
|------|--------|-------------|
| Stain Deconvolution | `stain_deconvolution.py` | Separates histological stains for improved image processing |
| Segmentation | `bin2cell-workflow.py` | Converts Visium HD bins into single-cell–like objects |
| Single-Cell Reference Generation | `generate_sc_reference.py` | Builds reference datasets for annotation |
| Cell Type Annotation | `tacco-annotate.py` | Annotates cell types using TACCO |
| Image Rotation / Alignment | `rotate-sample.py` | Adjusts orientation of spatial images |

---

## Coding Guidelines

No strict style guide was enforced during development.  
However, code formatting consistency is maintained using:

- **Black** (Python code formatter)

---

## Reproducibility & Environment Setup

We use **Conda** for dependency management.

### Requirements

- Install **Miniconda** or **Anaconda**

### Create the Environment

A predefined Conda environment file is included in the repository.  
You can recreate the environment by running:

```bash
make conda_env_create
