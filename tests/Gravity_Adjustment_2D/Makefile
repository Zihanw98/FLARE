CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile

.PHONY: test check
test: $(TARGET)
	python3 "$(CASE_DIR)/../../unit_test/gravity_checks.py" --case adjustment --run --executable "$(TARGET)" --parameters "$(CONFIG)"

check:
	python3 "$(CASE_DIR)/../../unit_test/gravity_checks.py" --case adjustment --parameters "$(CONFIG)"
