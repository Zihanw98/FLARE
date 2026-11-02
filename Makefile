# From any problem directory:
#   make -f /absolute/path/to/FLARE/Makefile -j4
#   make -f /absolute/path/to/FLARE/Makefile run
# Or from the FLARE root:
#   make PROBLEM_DIR=tests/sod_tube run
# The problem supplies user_problem.cpp (make_user_hooks) and user_config.h.

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

FLARE_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
FLARE_ROOT := $(patsubst %/,%,$(dir $(FLARE_MAKEFILE)))
PROBLEM_DIR ?= $(CURDIR)
override PROBLEM_DIR := $(abspath $(PROBLEM_DIR))
# Cache objects/dependencies separately; keep the executable next to the problem.
BUILD_DIR ?= $(PROBLEM_DIR)/.build
override BUILD_DIR := $(abspath $(BUILD_DIR))
TARGET ?= $(PROBLEM_DIR)/flare
override TARGET := $(abspath $(TARGET))
CONFIG ?= $(PROBLEM_DIR)/parameters.txt
override CONFIG := $(abspath $(CONFIG))

CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?=

# Collect every framework .cpp recursively, including main.cpp, but not old_version/tests.
find_cpp = $(wildcard $(1)*.cpp) $(foreach dir,$(wildcard $(1)*/),$(call find_cpp,$(dir)))
FRAMEWORK_SOURCES := $(sort $(call find_cpp,$(FLARE_ROOT)/src/))
# Flat, hidden object names encode source paths; framework/user names are separate.
framework_object = $(BUILD_DIR)/.flare-framework-$(subst /,__,$(patsubst $(FLARE_ROOT)/src/%.cpp,%,$(1))).o
user_object = $(BUILD_DIR)/.flare-user-$(subst /,__,$(patsubst %.cpp,%,$(1))).o
FRAMEWORK_OBJECTS := $(foreach source,$(FRAMEWORK_SOURCES),$(call framework_object,$(source)))

# Paths are relative to PROBLEM_DIR; override to select files or include subdirectories.
# Example: USER_SOURCES="user_problem.cpp helpers/initial_state.cpp"
USER_SOURCES ?= $(notdir $(wildcard $(PROBLEM_DIR)/*.cpp))
USER_OBJECTS := $(foreach source,$(USER_SOURCES),$(call user_object,$(source)))
OBJECTS := $(FRAMEWORK_OBJECTS) $(USER_OBJECTS)
ifneq ($(words $(OBJECTS)),$(words $(sort $(OBJECTS))))
$(error Source paths produce duplicate flat object names; rename the conflicting sources)
endif
DEPS := $(OBJECTS:.o=.d)

# The framework's config is a template, not a fallback include directory.
# Headers such as user/user_api.h remain accessible through the src include path.
FRAMEWORK_INCLUDE_DIRS := $(filter-out $(FLARE_ROOT)/src/user,$(sort $(patsubst %/,%,$(dir $(FRAMEWORK_SOURCES)))))
INCLUDE_FLAGS := -I"$(PROBLEM_DIR)" $(foreach dir,$(FRAMEWORK_INCLUDE_DIRS),-I"$(dir)")
override USER_CONFIG := $(PROBLEM_DIR)/user_config.h
override USER_CONFIG_FLAGS := -include "$(USER_CONFIG)"
# Also reject missing configs after setup, even when an old executable exists.
# Cleaning build products does not require a configuration header.
ifeq ($(wildcard $(USER_CONFIG)),)
ifneq ($(MAKECMDGOALS),clean)
$(error Missing problem configuration: $(USER_CONFIG). No configuration fallback is allowed)
endif
endif
BUILD_SETTINGS := $(BUILD_DIR)/.flare-build-settings

.PHONY: all run clean show-config check-problem FORCE
all: $(TARGET)

check-problem:
	@test -n "$(strip $(USER_SOURCES))" || { echo "No user .cpp found in $(PROBLEM_DIR). Set PROBLEM_DIR or USER_SOURCES." >&2; exit 1; }

# Rebuild on compiler/flag/source-selection changes.
# Leave the timestamp unchanged when the settings are identical (incremental builds).
$(BUILD_SETTINGS): FORCE
	@mkdir -p "$(@D)"
	@printf '%s\n' '$(CXX)' '$(CPPFLAGS)' '$(CXXFLAGS)' '$(LDFLAGS)' '$(LDLIBS)' '$(INCLUDE_FLAGS)' '$(USER_CONFIG_FLAGS)' '$(USER_CONFIG)' '$(OBJECTS)' > "$@.tmp"
	@if ! cmp -s "$@.tmp" "$@"; then mv "$@.tmp" "$@"; else $(RM) "$@.tmp"; fi

# Explicit source-to-object rules keep nested sources flat without copying them.
define compile_rule
$(1): $(2) $(USER_CONFIG) $(FLARE_MAKEFILE) $(BUILD_SETTINGS) | check-problem
	@mkdir -p "$$(@D)"
	$$(CXX) $$(INCLUDE_FLAGS) $$(CPPFLAGS) $$(CXXFLAGS) $$(USER_CONFIG_FLAGS) -MMD -MP -MF "$$(@:.o=.d)" -MT "$$@" -c "$$<" -o "$$@"
endef
$(foreach source,$(FRAMEWORK_SOURCES),$(eval $(call compile_rule,$(call framework_object,$(source)),$(source))))
$(foreach source,$(USER_SOURCES),$(eval $(call compile_rule,$(call user_object,$(source)),$(PROBLEM_DIR)/$(source))))

$(TARGET): $(OBJECTS) $(BUILD_SETTINGS)
	@mkdir -p "$(@D)"
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o "$@"

# Relative output paths in parameters.txt are relative to the user's problem directory.
run: $(TARGET)
	cd "$(PROBLEM_DIR)" && "$(TARGET)" -par "$(CONFIG)"

show-config:
	@printf '%s\n' 'FLARE_ROOT=$(FLARE_ROOT)' 'PROBLEM_DIR=$(PROBLEM_DIR)' 'BUILD_DIR=$(BUILD_DIR)' 'USER_CONFIG=$(USER_CONFIG)' 'USER_SOURCES=$(USER_SOURCES)' 'TARGET=$(TARGET)' 'CONFIG=$(CONFIG)'

# Remove only this build's products; preserve simulation CSV files and user sources.
clean:
	@$(RM) $(OBJECTS) $(DEPS) "$(TARGET)" "$(BUILD_SETTINGS)" "$(BUILD_SETTINGS).tmp"

-include $(DEPS)
