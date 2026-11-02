FLARE_TRANSLATION_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(FLARE_TRANSLATION_DIR)
include $(FLARE_TRANSLATION_DIR)/../../Makefile

.PHONY: test check animate
# Run a simulation, validate translation/conservation, and save comparison.png.
test: $(TARGET)
	python3 "$(FLARE_TRANSLATION_DIR)/check.py" --run --executable "$(TARGET)" --parameters "$(CONFIG)"

# Analyze only files named by the existing run.log; never mix old snapshots in.
check:
	python3 "$(FLARE_TRANSLATION_DIR)/check.py" --parameters "$(CONFIG)"

animate:
	python3 "$(FLARE_TRANSLATION_DIR)/check.py" --parameters "$(CONFIG)" --animate
