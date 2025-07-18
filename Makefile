# Makefile for neutralb1 project
# This makefile compiles C++ source files into the build directory structure

# Build configuration (debug or release)
BUILD_TYPE ?= release

# Project environment file
ENV_FILE := config/.env
SETUP_SCRIPT := config/setup_gluex.sh
-include $(ENV_FILE)

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
    BASE_BUILD_CXXFLAGS := $(BASE_CXXFLAGS) -g -O0 -DDEBUG
else ifeq ($(BUILD_TYPE),release)
    BASE_BUILD_CXXFLAGS := $(BASE_CXXFLAGS) -O2 -DNDEBUG
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
ROOT_CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOT_LIBS = $(shell $(ROOTSYS)/bin/root-config --libs)

# Library paths from your setup script
LIBRARY_DIRS := -L$(ROOTSYS)/lib \
                -L$(AMPTOOLS)/lib \
                -L$(HALLD_SIM_HOME)/Linux_Alma9-x86_64-gcc11.5.0/lib

# Libraries to link against
LIBS = $(ROOT_LIBS) -lAmpTools -lAMPTOOLS_AMPS

# Combined flags
CXXFLAGS = $(BASE_BUILD_CXXFLAGS) $(ROOT_CFLAGS) $(INCLUDE_DIRS)
LDFLAGS = $(LIBRARY_DIRS) $(LIBS)

# Default target
.PHONY: all
all: $(ENV_FILE) directories $(BATCH_EXECUTABLES) $(ALL_OBJECTS)

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

# Rule to generate the environment file from the setup script
$(ENV_FILE): $(SETUP_SCRIPT)
	@echo "Generating environment file from setup script..."
	@mkdir -p config
	@bash -c 'source $(SETUP_SCRIPT) && env | grep -E "^(ROOTSYS|AMPTOOLS|HALLD_|FSROOT|ROOT_|PATH|LD_LIBRARY_PATH|LIBRARY_PATH|CPLUS_INCLUDE_PATH|DYLD_LIBRARY_PATH)" | sort > $(ENV_FILE)'
	@echo "PROJECT_ROOT=$(PROJECT_ROOT)" >> $(ENV_FILE)
	@echo "Environment file generated at $(ENV_FILE)"

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
	@echo "Cleaning all build directories and the environment file..."
	rm -rf build/*
	rm -f $(ENV_FILE)

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

.PHONY: clean-env
clean-env:
	@echo "Cleaning environment file..."
	rm -f $(ENV_FILE)

.PHONY: update-env
update-env:
	@echo "Forcing regeneration of environment file..."
	@rm -f $(ENV_FILE)
	@$(MAKE) $(ENV_FILE)


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
	@echo "  clean            - Remove all build files and environment file"
	@echo "  clean-debug      - Remove only debug build files"
	@echo "  clean-release    - Remove only release build files"
	@echo "  clean-objects    - Remove only object files"
	@echo "  clean-bins       - Remove only executables"
	@echo "  clean-env        - Remove environment file"	
	@echo "  update-env       - Force regeneration of environment file"
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
	@echo ""
	@echo "Environment:"
	@echo "  The .env file is automatically generated from $(SETUP_SCRIPT)"
	@echo "  Use 'make update-env' to force regeneration"

# Dependency tracking (optional, for automatic rebuilds when headers change)
# Only include dependencies if we're not cleaning and environment is loaded
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),clean-all)
ifneq ($(MAKECMDGOALS),clean-env)
ifneq ($(MAKECMDGOALS),update-env)
ifdef ROOTSYS
-include $(ALL_OBJECTS:.o=.d)
endif
endif
endif
endif
endif

$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT $(@:.d=.o) $< > $@

.PHONY: deps
deps: $(ALL_OBJECTS:.o=.d)
