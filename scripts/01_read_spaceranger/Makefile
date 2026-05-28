.PHONY: read

ENV_NAME ?= spatialenv_publication
SAMPLES ?=
BASE_OUTPUT ?= results
CONFIG_DIR ?= config

# repo root = parent of scripts/
REPO_ROOT := $(abspath $(CURDIR)/../..)

PYTHON ?= conda run -n $(ENV_NAME) python

read:
	@mkdir -p $(REPO_ROOT)/$(BASE_OUTPUT)/zarr_before_qc
	@for s in $(SAMPLES); do \
		echo "==> Reading SpaceRanger for $$s"; \
		cd $(REPO_ROOT) && \
		$(PYTHON) scripts/01_read_spaceranger/read_spaceRangerv4.py \
			--config $(CONFIG_DIR)/$$s.yaml \
			--output-dir $(BASE_OUTPUT); \
	done
