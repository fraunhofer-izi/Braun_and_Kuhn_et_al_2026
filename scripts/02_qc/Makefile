.PHONY: qc

ENV_NAME ?= spatialenv_publication
SAMPLES ?=
BASE_OUTPUT ?= results
CONFIG_DIR ?= config

# repo root = two levels up from scripts/02_qc
REPO_ROOT := $(abspath $(CURDIR)/../..)

PYTHON ?= conda run -n $(ENV_NAME) python

qc:
	@mkdir -p $(REPO_ROOT)/$(BASE_OUTPUT)/zarr_after_qc
	@for s in $(SAMPLES); do \
		echo "==> QC for $$s"; \
		cd $(REPO_ROOT) && \
		$(PYTHON) scripts/02_qc/spatial_data_qc.py \
			--config $(CONFIG_DIR)/$$s.yaml \
			--zarr-path $(BASE_OUTPUT)/zarr_before_qc/$$s.zarr \
			--output-dir $(BASE_OUTPUT); \
	done
