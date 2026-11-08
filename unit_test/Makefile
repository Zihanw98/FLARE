# Independent executables; no simulation main or user_problem.cpp is linked.
.DEFAULT_GOAL := test
.DELETE_ON_ERROR:
THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
UNIT_DIR := $(patsubst %/,%,$(dir $(THIS_MAKEFILE)))
ROOT := $(abspath $(UNIT_DIR)/..)
CXX ?= g++
AR ?= ar
CXXFLAGS ?= -O1 -g -std=c++17 -Wall -Wextra -pedantic
DIMENSIONS ?= 1 2 3
NDIM ?= 1
AVAILABLE_TESTS := state mesh facefield godunov integrator rk2 source wave_timestep output_boundary limiter units initialization performance
TESTS ?= $(AVAILABLE_TESTS)
ifneq ($(filter-out $(AVAILABLE_TESTS),$(TESTS)),)
$(error Unknown TESTS: $(filter-out $(AVAILABLE_TESTS),$(TESTS)))
endif
ifneq ($(filter-out 1 2 3,$(DIMENSIONS) $(NDIM)),)
$(error DIMENSIONS and NDIM must be 1, 2 or 3)
endif

BUILD_DIR := $(UNIT_DIR)/build/ndim$(NDIM)
find_cpp = $(wildcard $(1)*.cpp) $(foreach dir,$(wildcard $(1)*/),$(call find_cpp,$(dir)))
SOURCES := $(filter-out $(ROOT)/src/main.cpp,$(sort $(call find_cpp,$(ROOT)/src/)))
OBJECTS := $(patsubst $(ROOT)/src/%.cpp,$(BUILD_DIR)/framework/%.o,$(SOURCES))
INCLUDES := $(foreach dir,$(sort $(dir $(SOURCES))) $(ROOT)/src $(ROOT)/src/user,-I"$(dir)")
# Old tests use assert: always enable assertions, even if CXXFLAGS contains -DNDEBUG.
FLAGS := $(INCLUDES) $(CPPFLAGS) $(CXXFLAGS) -DHD_NDIM=$(NDIM) -UNDEBUG
LIBRARY := $(BUILD_DIR)/libflare.a
BINARIES := $(addprefix $(BUILD_DIR)/test_,$(TESTS))
SETTINGS := $(BUILD_DIR)/.build-settings
vpath test_%.cpp $(UNIT_DIR) $(ROOT)/tests

.PHONY: test check FORCE
test:
	+@set -eu; for dimension in $(DIMENSIONS); do \
	    $(MAKE) --no-print-directory -f "$(THIS_MAKEFILE)" check NDIM=$$dimension; \
	done

# Each run gets disposable output/config files; no Sod snapshots are touched.
check: $(BINARIES)
	@set -eu; test_dir=$$(mktemp -d); trap 'rm -rf "$$test_dir"' EXIT; \
	for name in $(TESTS); do \
	    printf '\n[NDIM=$(NDIM)] %s\n' "$$name"; \
	    "$(BUILD_DIR)/test_$$name" "$$test_dir/$$name"; \
	done

$(SETTINGS): FORCE
	@mkdir -p "$(@D)"
	@printf '%s\n' '$(CXX)' '$(AR)' '$(FLAGS)' '$(LDFLAGS)' '$(LDLIBS)' '$(OBJECTS)' > "$@.tmp"
	@if ! cmp -s "$@.tmp" "$@"; then mv "$@.tmp" "$@"; else $(RM) "$@.tmp"; fi

$(OBJECTS): $(SETTINGS) $(THIS_MAKEFILE)
$(BUILD_DIR)/framework/%.o: $(ROOT)/src/%.cpp
	@mkdir -p "$(@D)"
	$(CXX) $(FLAGS) -MMD -MP -c "$<" -o "$@"

$(LIBRARY): $(OBJECTS) $(SETTINGS)
	@$(RM) "$@"
	$(AR) rcs "$@" $(OBJECTS)

$(BUILD_DIR)/test_%: test_%.cpp $(LIBRARY) $(SETTINGS) $(THIS_MAKEFILE)
	$(CXX) $(FLAGS) -MMD -MP -MF "$@.d" -MT "$@" "$<" "$(LIBRARY)" $(LDFLAGS) $(LDLIBS) -o "$@"

-include $(OBJECTS:.o=.d) $(BINARIES:=.d)
