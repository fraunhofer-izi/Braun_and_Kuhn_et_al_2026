.PHONY: download_data link_duodenum

# Wo liegt dein echtes SpaceRanger outs-Verzeichnis (AUßERHALB des Repos)?
# beim Aufruf als Variable setzen:
#   make 00_download_data SRC_OUTS=/abs/pfad/.../outs
SRC_OUTS ?= 

# Ziel im Repo:
DEST_DIR := data/spaceranger/Duodenum_run3
DEST_OUTS := $(DEST_DIR)/outs

download_data: link_duodenum

link_duodenum:
	@if [ -z "$(SRC_OUTS)" ]; then \
		echo "ERROR: Bitte SRC_OUTS setzen (absoluter Pfad zum echten outs)."; \
		echo "Beispiel: make 00_download_data SRC_OUTS=/mnt/.../Duodenum_run3/outs"; \
		exit 1; \
	fi
	@mkdir -p $(DEST_DIR)
	@if [ -e "$(DEST_OUTS)" ] || [ -L "$(DEST_OUTS)" ]; then \
		echo "INFO: Entferne vorhandenes $(DEST_OUTS)"; \
		rm -rf "$(DEST_OUTS)"; \
	fi
	@ln -s "$(SRC_OUTS)" "$(DEST_OUTS)"
	@echo "OK: Symlink erstellt:"
	@ls -la "$(DEST_OUTS)"
