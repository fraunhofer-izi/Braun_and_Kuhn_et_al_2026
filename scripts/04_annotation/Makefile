.PHONY: anno

ENV_NAME ?= spatialenv_publication
SAMPLES ?=
BASE_OUTPUT ?= results
CONFIG_DIR ?= config

REPO_ROOT := $(abspath $(CURDIR)/../..)
PYTHON ?= conda run -n $(ENV_NAME) python

anno:
	@mkdir -p $(REPO_ROOT)/$(BASE_OUTPUT)/zarr_after_anno
	@for s in $(SAMPLES); do \
		echo "==> Annotation $$s"; \
		cd $(REPO_ROOT) && \
		$(PYTHON) scripts/04_annotation/annotation_with_aucell.py \
			--config $(CONFIG_DIR)/$$s.yaml \
			--zarr-path $(BASE_OUTPUT)/zarr_after_cluster/$$s.zarr \
			--output-dir $(BASE_OUTPUT); \
	done
