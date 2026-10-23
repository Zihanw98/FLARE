CXX = g++
CPPFLAGS += -Isrc
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic
BUILD_DIR = build
TARGET = $(BUILD_DIR)/hdbasic
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPS = $(OBJECTS:.o=.d)
CONFIG ?= parameters.txt

.PHONY: all run clean
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: src/%.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR):
	mkdir -p $@

# Snapshots are written relative to the process working directory.
run: $(TARGET)
	mkdir -p output
	cd output && "$(abspath $(TARGET))" "$(abspath $(CONFIG))"

# Keep simulation snapshots and figures.
clean:
	$(RM) $(OBJECTS) $(DEPS) $(TARGET)

-include $(DEPS)
