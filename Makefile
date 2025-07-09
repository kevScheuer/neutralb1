# Makefile for neutralb1 project
# This makefile compiles C++ source files into the build directory structure

# Build configuration (debug or release)
BUILD_TYPE ?= release

# Project directories
PROJECT_ROOT := $(shell pwd)
SRC_DIR := src
BUILD_DIR := build/$(BUILD_TYPE)
BIN_DIR := $(BUILD_DIR)/bin
LIB_DIR := $(BUILD_DIR)/lib
OBJ_DIR := $(BUILD_DIR)/obj

# Source files
BATCH_SOURCES := $(wildcard $(SRC_DIR)/batch/*.cc)
ANALYSIS_SOURCES := $(wildcard $(SRC_DIR)/analysis/*.cc)
SELECTION_SOURCES := $(wildcard $(SRC_DIR)/selection/*.cc)
UTILS_SOURCES := $(wildcard $(SRC_DIR)/utils/*.cc)

# Object files (preserving directory structure)
BATCH_OBJECTS := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(BATCH_SOURCES))
ANALYSIS_OBJECTS := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(ANALYSIS_SOURCES))
SELECTION_OBJECTS := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SELECTION_SOURCES))
UTILS_OBJECTS := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(UTILS_SOURCES))

ALL_OBJECTS := $(BATCH_OBJECTS) $(ANALYSIS_OBJECTS) $(SELECTION_OBJECTS) $(UTILS_OBJECTS)

# Executables (batch tools become binaries)
BATCH_EXECUTABLES := $(patsubst $(SRC_DIR)/batch/%.cc,$(BIN_DIR)/%,$(BATCH_SOURCES))

# Compiler settings
CXX := g++
BASE_CXXFLAGS := -std=c++17 -Wall -Wextra -fPIC

# Build-specific flags
ifeq ($(BUILD_TYPE),debug)
    CXXFLAGS := $(BASE_CXXFLAGS) -g -O0 -DDEBUG
else ifeq ($(BUILD_TYPE),release)
    CXXFLAGS := $(BASE_CXXFLAGS) -O2 -DNDEBUG
else
    $(error Invalid BUILD_TYPE: $(BUILD_TYPE). Use 'debug' or 'release')
endif

# Include paths from your setup script
INCLUDE_DIRS := -I$(PROJECT_ROOT)/include \
                -I$(ROOTSYS)/include \
                -I$(AMPTOOLS) \
                -I$(HALLD_SIM_HOME)/src/libraries \
                -I$(HALLD_SIM_HOME)/Linux_Alma9-x86_64-gcc11.5.0/include

# ROOT specific flags
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)

# Library paths from your setup script
LIBRARY_DIRS := -L$(ROOTSYS)/lib \
                -L$(AMPTOOLS)/lib \
                -L$(HALLD_SIM_HOME)/Linux_Alma9-x86_64-gcc11.5.0/lib

# Libraries to link against
LIBS := $(ROOT_LIBS) -lAmpTools

# Combined flags
CXXFLAGS += $(ROOT_CFLAGS) $(INCLUDE_DIRS)
LDFLAGS := $(LIBRARY_DIRS) $(LIBS)

# Default target
.PHONY: all
all: directories $(BATCH_EXECUTABLES) $(ALL_OBJECTS)

# Debug and release targets
.PHONY: debug release
debug:
	$(MAKE) BUILD_TYPE=debug all

release:
	$(MAKE) BUILD_TYPE=release all

# Create build directories
.PHONY: directories
directories:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(OBJ_DIR)/batch
	@mkdir -p $(OBJ_DIR)/analysis
	@mkdir -p $(OBJ_DIR)/selection
	@mkdir -p $(OBJ_DIR)/utils

# Rule to compile object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@echo "Compiling $<..."
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to create executables from batch sources
$(BIN_DIR)/%: $(OBJ_DIR)/batch/%.o $(UTILS_OBJECTS) $(ANALYSIS_OBJECTS) $(SELECTION_OBJECTS)
	@echo "Linking $@..."
	$(CXX) $^ -o $@ $(LDFLAGS)

# Individual executable targets for convenience
.PHONY: extract_fit_results angle_plotter extract_cov_matrix extract_corr_matrix extract_bin_info project_moments
extract_fit_results: $(BIN_DIR)/extract_fit_results
angle_plotter: $(BIN_DIR)/angle_plotter
extract_cov_matrix: $(BIN_DIR)/extract_cov_matrix
extract_corr_matrix: $(BIN_DIR)/extract_corr_matrix
extract_bin_info: $(BIN_DIR)/extract_bin_info
project_moments: $(BIN_DIR)/project_moments

# Create a shared library from utility and analysis objects
$(LIB_DIR)/libneutralb1.so: $(UTILS_OBJECTS) $(ANALYSIS_OBJECTS) $(SELECTION_OBJECTS)
	@echo "Creating shared library $@..."
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

.PHONY: library
library: $(LIB_DIR)/libneutralb1.so

# Clean targets
.PHONY: clean
clean:
	@echo "Cleaning all build directories..."
	rm -rf build/*

.PHONY: clean-debug
clean-debug:
	@echo "Cleaning debug build directory..."
	rm -rf build/debug/*

.PHONY: clean-release
clean-release:
	@echo "Cleaning release build directory..."
	rm -rf build/release/*

.PHONY: clean-objects
clean-objects:
	@echo "Cleaning object files..."
	rm -rf $(OBJ_DIR)/*

.PHONY: clean-bins
clean-bins:
	@echo "Cleaning executables..."
	rm -rf $(BIN_DIR)/*

# Setup environment (sources your setup script)
.PHONY: setup
setup:
	@echo "To setup environment, source the setup script:"
	@echo "source config/setup_gluex.sh"

# Install targets (copies binaries to a system location if needed)
INSTALL_PREFIX ?= /usr/local
.PHONY: install
install: all
	@echo "Installing to $(INSTALL_PREFIX)..."
	mkdir -p $(INSTALL_PREFIX)/bin
	cp $(BIN_DIR)/* $(INSTALL_PREFIX)/bin/
	mkdir -p $(INSTALL_PREFIX)/lib
	cp $(LIB_DIR)/* $(INSTALL_PREFIX)/lib/ 2>/dev/null || true

# Debug information
.PHONY: info
info:
	@echo "Project Root: $(PROJECT_ROOT)"
	@echo "Source files found:"
	@echo "  Batch: $(BATCH_SOURCES)"
	@echo "  Analysis: $(ANALYSIS_SOURCES)"
	@echo "  Selection: $(SELECTION_SOURCES)"
	@echo "  Utils: $(UTILS_SOURCES)"
	@echo ""
	@echo "Executables will be built:"
	@echo "  $(BATCH_EXECUTABLES)"
	@echo ""
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Libraries: $(LDFLAGS)"

# Help target
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all              - Build all executables and object files (current: $(BUILD_TYPE))"
	@echo "  debug            - Build debug version with -g -O0"
	@echo "  release          - Build release version with -O2"
	@echo "  clean            - Remove all build files"
	@echo "  clean-debug      - Remove only debug build files"
	@echo "  clean-release    - Remove only release build files"
	@echo "  clean-objects    - Remove only object files"
	@echo "  clean-bins       - Remove only executables"
	@echo "  library          - Build shared library from analysis/utils code"
	@echo "  setup            - Show setup instructions"
	@echo "  install          - Install binaries to $(INSTALL_PREFIX)"
	@echo "  info             - Show project information"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Individual executables:"
	@echo "  extract_fit_results"
	@echo "  angle_plotter"
	@echo "  extract_cov_matrix"
	@echo "  extract_corr_matrix"
	@echo "  extract_bin_info"
	@echo "  project_moments"

# Dependency tracking (optional, for automatic rebuilds when headers change)
-include $(ALL_OBJECTS:.o=.d)

$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT $(@:.d=.o) $< > $@

.PHONY: deps
deps: $(ALL_OBJECTS:.o=.d)
