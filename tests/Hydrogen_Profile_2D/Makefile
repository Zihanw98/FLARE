CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
PYTHON ?= python3
include $(CASE_DIR)/../../Makefile

.PHONY: inputs test animate
# Regenerate only on explicit request, so normal runs preserve manual profile edits.
inputs:
	$(PYTHON) "$(CASE_DIR)/generate_profile.py"

test: run

animate:
	$(PYTHON) "$(CASE_DIR)/animate.py" --parameters "$(CONFIG)"
