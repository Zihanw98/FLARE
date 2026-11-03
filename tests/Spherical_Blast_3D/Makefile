CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile

.PHONY: test check animate render
test: $(TARGET)
	python3 "$(CASE_DIR)/check.py" --run --executable "$(TARGET)" --parameters "$(CONFIG)"

check:
	python3 "$(CASE_DIR)/check.py" --parameters "$(CONFIG)"

# The manifest in output/run.log excludes snapshots left by earlier runs.
animate:
	python3 "$(CASE_DIR)/../../animate3d.py" "$(CASE_DIR)/output" mesh_blast3d --field p

render:
	python3 "$(CASE_DIR)/../../render3d.py" "$(CASE_DIR)/output" mesh_blast3d --obj-sequence
