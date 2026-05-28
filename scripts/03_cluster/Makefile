.PHONY: cluster

ENV_NAME ?= spatialenv_publication
SAMPLES ?=
BASE_OUTPUT ?= results
CONFIG_DIR ?= config

REPO_ROOT := $(abspath $(CURDIR)/../..)
PYTHON ?= conda run -n $(ENV_NAME) python

cluster:
	@mkdir -p $(REPO_ROOT)/$(BASE_OUTPUT)/zarr_after_cluster
	@for s in $(SAMPLES); do \
		echo "==> Clustering $$s"; \
		cd $(REPO_ROOT) && \
		$(PYTHON) scripts/03_cluster/norm_and_cluster.py \
			--config $(CONFIG_DIR)/$$s.yaml \
			--zarr-path $(BASE_OUTPUT)/zarr_after_qc/$$s.zarr \
			--output-dir $(BASE_OUTPUT); \
	done
