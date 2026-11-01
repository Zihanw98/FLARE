# make / make run / make test, from this directory or using make -C tests/sod_tube.
SOD_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(SOD_DIR)
include $(SOD_DIR)/../../Makefile

.PHONY: test animate
test: $(TARGET)
	python3 "$(SOD_DIR)/check.py" --executable "$(TARGET)" --parameters "$(CONFIG)"

# Read existing snapshots; `make test` generates CSVs and validates the run.
animate:
	python3 "$(SOD_DIR)/animate.py" $(ANIMATION_ARGS)
