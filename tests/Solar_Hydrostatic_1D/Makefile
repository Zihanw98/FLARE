CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile

.PHONY: test check
test: $(TARGET)
	python3 "$(CASE_DIR)/check.py" --run --parameters "$(CONFIG)" --executable "$(TARGET)"

check:
	python3 "$(CASE_DIR)/check.py" --parameters "$(CONFIG)"
