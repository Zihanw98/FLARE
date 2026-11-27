CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
USER_SOURCES := ../../unit_test/snec_thermal_bomb.cpp
include $(CASE_DIR)/../../Makefile

.PHONY: compare analyze
compare: $(TARGET)
	python3 $(CASE_DIR)/../../unit_test/compare_snec_thermal_bomb.py run

analyze:
	python3 $(CASE_DIR)/../../unit_test/compare_snec_thermal_bomb.py analyze
