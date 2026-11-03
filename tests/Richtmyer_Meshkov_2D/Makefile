CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile
.PHONY: test check animate
test: $(TARGET)
	python3 "$(CASE_DIR)/../hydro_suite.py" rm --run --executable "$(TARGET)" --parameters "$(CONFIG)"
check:
	python3 "$(CASE_DIR)/../hydro_suite.py" rm --parameters "$(CONFIG)"
animate:
	python3 "$(CASE_DIR)/../hydro_suite.py" rm --parameters "$(CONFIG)" --animate
