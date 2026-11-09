CASE_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PROBLEM_DIR ?= $(CASE_DIR)
include $(CASE_DIR)/../../Makefile

CHECKER := $(BUILD_DIR)/.flare-check-residual
CHECK_OBJECTS := $(filter-out $(call framework_object,$(FLARE_ROOT)/src/main.cpp),$(OBJECTS))
$(CHECKER): $(CASE_DIR)/verification/check_residual.cpp $(CHECK_OBJECTS) $(BUILD_SETTINGS)
	$(CXX) $(INCLUDE_FLAGS) $(CPPFLAGS) $(CXXFLAGS) $(USER_CONFIG_FLAGS) $< $(CHECK_OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

.PHONY: check-residual
check-residual: $(CHECKER)
	cd "$(CASE_DIR)" && "$(CHECKER)" "$(CONFIG)"
