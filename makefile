PROJECT_ROOT := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
EXAMPLES_DIR := $(PROJECT_ROOT)examples
SAGE_FILES   := $(shell find $(PROJECT_ROOT) -type f -name '*.sage')
BACKUP_DIR   := $(PROJECT_ROOT)backups

.PHONY: all write-root restore clean uninstall

all: write-root

write-root:
	@mkdir -p "$(BACKUP_DIR)"
	@echo "Injecting PROJECT_ROOT='$(PROJECT_ROOT)' into .sage files"
	@set -e; \
	INJ="PROJECT_ROOT = '$(PROJECT_ROOT)'"; \
	for f in $(SAGE_FILES); do \
		base=$$(basename "$$f"); \
		echo "  -> $$f"; \
		cp "$$f" "$(BACKUP_DIR)/$$base.bak"; \
		if head -n1 "$$f" | grep -q '^#!'; then \
			{ head -n1 "$$f"; echo "$$INJ"; echo; tail -n +2 "$$f" | sed '/^PROJECT_ROOT *=/d'; } > "$$f.tmp"; \
		else \
			{ echo "$$INJ"; echo; sed '/^PROJECT_ROOT *=/d' "$$f"; } > "$$f.tmp"; \
		fi; \
		mv "$$f.tmp" "$$f"; \
	done
	@echo "Done."


restore:
	@echo "Restoring .sage files from backups/"
	@set -e; \
	for b in $(shell find "$(BACKUP_DIR)" -type f -name '*.bak'); do \
		base=$$(basename "$$b" .bak); \
		target=$$(find "$(PROJECT_ROOT)" -type f -name "$$base" | head -n1); \
		if [ -n "$$target" ]; then \
			echo "  <- $$target"; \
			cp "$$b" "$$target"; \
		fi; \
	done
	@echo "Done."

clean:
	@echo "Removing backups/"
	@rm -rf "$(BACKUP_DIR)"

uninstall: restore clean
	@echo "Uninstalled injected roots and removed backups."
