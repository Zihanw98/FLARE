CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile

.PHONY: inputs test
inputs:
	python3 "$(CASE_DIR)/generate_inputs.py"

test: $(TARGET)
	cd "$(CASE_DIR)" && "$(TARGET)" -par parameters.txt > run.log
	python3 "$(CASE_DIR)/check.py"
