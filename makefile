PROJECT_ROOT := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ROOT_FILE_MAIN := $(PROJECT_ROOT)PROJECT_ROOT
ROOT_FILE_EXAMPLES := $(PROJECT_ROOT)examples/PROJECT_ROOT

EXAMPLES_DIR := $(PROJECT_ROOT)examples


.PHONY: all write-root uninstall 

all: write-root 
	

write-root:
	@mkdir -p "$(EXAMPLES_DIR)"
	@echo "$(PROJECT_ROOT)" > "$(ROOT_FILE_MAIN)"
	@echo "$(PROJECT_ROOT)" > "$(ROOT_FILE_EXAMPLES)"
	@echo "Set PROJECT_ROOT to: $(PROJECT_ROOT)"

uninstall:
	@rm -f "$(ROOT_FILE_MAIN)" "$(ROOT_FILE_EXAMPLES)"
	@echo "Removed Root Files"
	@if [ -d "$(EXAMPLES_DIR)/metadata" ]; then \
		rm -rf "$(EXAMPLES_DIR)/metadata" && \
		echo "Removed $(EXAMPLES_DIR)/metadata"; \
		fi; 