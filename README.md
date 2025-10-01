# Introduction

 :bangbang:**This repository must be on a JLab ifarm node**:bangbang:

Within this repository lies all scripts needed to perform a full partial wave analysis of the neutral $\omega\pi^0$ channel from start to finish, which is dominated by a neutral $b_1$ resonance for which this repo is named after. 

# Setup

This project uses `uv` with `scikit-build-core` and CMake to build both Python packages and C++ executables. The build system supports both debug and release configurations.

## Prerequisites

1. **JLab ifarm node** - This repository must be on a JLab ifarm node
2. **uv** - Python package manager ([installation instructions](https://docs.astral.sh/uv/getting-started/installation/))
3. **GlueX environment** - Source the environment first: `source config/setup_gluex.sh`

## Build Options

Use the provided build script for convenience:

```bash
# Quick builds for c++ binaries - outputs to build/debug/ or build/release/
./build_with_env.sh debug    # Debug build → build/debug/bin/
./build_with_env.sh release  # Release build → build/release/bin/

# Python package builds - outputs to build/python/
./build_with_env.sh dev Debug    # Development debug build
./build_with_env.sh dev Release  # Development release build

# Wheel builds - outputs to build/python/ and dist/
./build_with_env.sh wheel Debug    # Debug wheel
./build_with_env.sh wheel Release  # Release wheel (default)
```

## Build Outputs

- **`debug`/`release`**: Executables in `build/debug/bin/` or `build/release/bin/`
- **`dev`**: Python package installed, C++ executables accessible via Python environment
- **`wheel`**: Self-contained Python wheel with embedded C++ executables

## Manual Build

If you prefer to build manually:

```bash
# Source the environment first
source config/setup_gluex.sh

# Set build type (Debug or Release)
export CMAKE_BUILD_TYPE=Release

# Install in development mode
uv pip install -e .

# Or build a wheel
uv build
```

# Documentation
Documentation is currently compiled by running `make html` within the [`docs` folder](./docs/).

TODO: Describe where to find documentation, and maybe a short description of how to make it. Might want to use sphinx + doxygen [as described here](https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/), but need to look into other suggestions for easiest solution to making dual c++/python docs