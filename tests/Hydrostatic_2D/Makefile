CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile

.PHONY: test check diagnose
# Strict equilibrium preservation: expected to fail with the present unbalanced scheme.
test: $(TARGET)
	python3 "$(CASE_DIR)/../gravity_checks.py" --case hydrostatic --run --executable "$(TARGET)" --parameters "$(CONFIG)"

check:
	python3 "$(CASE_DIR)/../gravity_checks.py" --case hydrostatic --parameters "$(CONFIG)"

# Run all invariant checks and report equilibrium drift without failing on that drift.
diagnose: $(TARGET)
	python3 "$(CASE_DIR)/../gravity_checks.py" --case hydrostatic --diagnostic-only --run --executable "$(TARGET)" --parameters "$(CONFIG)"
