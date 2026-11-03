CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile
.PHONY: test check animate converge
test: $(TARGET)
	python3 "$(CASE_DIR)/../hydro_suite.py" acoustic --run --executable "$(TARGET)" --parameters "$(CONFIG)"
check:
	python3 "$(CASE_DIR)/../hydro_suite.py" acoustic --parameters "$(CONFIG)"
animate:
	python3 "$(CASE_DIR)/../hydro_suite.py" acoustic --parameters "$(CONFIG)" --animate
converge: $(TARGET)
	python3 "$(CASE_DIR)/../hydro_suite.py" acoustic --convergence --executable "$(TARGET)" --parameters "$(CONFIG)"
