CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile

.PHONY: test check animate
test: $(TARGET)
	python3 "$(CASE_DIR)/check.py" --run --executable "$(TARGET)" --parameters "$(CONFIG)"

check:
	python3 "$(CASE_DIR)/check.py" --parameters "$(CONFIG)"

animate:
	python3 "$(CASE_DIR)/check.py" --parameters "$(CONFIG)" --animate
