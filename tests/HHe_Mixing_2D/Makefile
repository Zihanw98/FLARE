CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
PYTHON ?= python3
include $(CASE_DIR)/../../Makefile
.PHONY: test check animate smoke
test: $(TARGET)
	$(PYTHON) "$(CASE_DIR)/check.py" --run --executable "$(TARGET)" --parameters "$(CONFIG)"
check:
	$(PYTHON) "$(CASE_DIR)/check.py" --parameters "$(CONFIG)"
animate:
	$(PYTHON) "$(CASE_DIR)/animate.py" --parameters "$(CONFIG)"
smoke: $(TARGET)
	$(PYTHON) "$(CASE_DIR)/check.py" --run --smoke --executable "$(TARGET)" --parameters "$(CONFIG)"
