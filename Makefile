.PHONY: conda_env_create 00_download_data 01_read 02_qc 03_cluster 04_anno all

# User-overridable settings
SAMPLES ?= Skin Duodenum_run3
BASE_OUTPUT ?= results
ENV_FILE ?= environment.yml
ENV_NAME ?= spatialenv_publication

conda_env_create: $(ENV_FILE)
	mamba env create -f $(ENV_FILE) -n $(ENV_NAME) || mamba env update -f $(ENV_FILE) -n $(ENV_NAME)
	$(MAKE) install_local_package

install_local_package:
	conda run -n $(ENV_NAME) python -m pip install -e .

00_download_data: scripts/00_download_data/Makefile
	$(MAKE) -C scripts/00_download_data download_data

01_read: scripts/01_read_spaceranger/Makefile
	$(MAKE) -C scripts/01_read_spaceranger read ENV_NAME=$(ENV_NAME) SAMPLES="$(SAMPLES)" BASE_OUTPUT=$(BASE_OUTPUT)

02_qc: scripts/02_qc/Makefile
	$(MAKE) -C scripts/02_qc qc ENV_NAME=$(ENV_NAME) SAMPLES="$(SAMPLES)" BASE_OUTPUT=$(BASE_OUTPUT)

03_cluster: scripts/03_cluster/Makefile
	$(MAKE) -C scripts/03_cluster cluster ENV_NAME=$(ENV_NAME) SAMPLES="$(SAMPLES)" BASE_OUTPUT=$(BASE_OUTPUT)

04_anno: scripts/04_annotation/Makefile
	$(MAKE) -C scripts/04_annotation anno ENV_NAME=$(ENV_NAME) SAMPLES="$(SAMPLES)" BASE_OUTPUT=$(BASE_OUTPUT)

all: 00_download_data 01_read 02_qc 03_cluster 04_anno
