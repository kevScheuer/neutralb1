#!/bin/bash
# build_with_env.sh - Helper script to build with proper environment

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

# Add local bin to PATH to find uv
export PATH="$HOME/.local/bin:$PATH"
source $PROJECT_ROOT/.venv/bin/activate

# Source the GlueX environment and create an env file of the necessary variables
if [ -f "$PROJECT_ROOT/config/setup_gluex.sh" ]; then
    ENV_FILE="$PROJECT_ROOT/config/.env"
    echo "Sourcing GlueX environment..."
    source "$PROJECT_ROOT/config/setup_gluex.sh" "$PROJECT_ROOT/config/version.xml"
    env | grep -E "^(ROOTSYS|AMPTOOLS|HALLD_|FSROOT|ROOT_|PATH|LD_LIBRARY_PATH|LIBRARY_PATH|CPLUS_INCLUDE_PATH|DYLD_LIBRARY_PATH)" | sort > "$ENV_FILE"
else
    echo "Warning: GlueX setup script not found. Build may fail if dependencies are not available."
fi

# Check if we're in development mode or building a wheel
if [ "$1" = "wheel" ]; then
    BUILD_TYPE="${2:-Release}"
    echo "Building wheel with uv and scikit-build-core (${BUILD_TYPE})..."
    cd "$PROJECT_ROOT"
    
    # Set build type for scikit-build-core
    export CMAKE_BUILD_TYPE="$BUILD_TYPE"
    uv build
    
elif [ "$1" = "dev" ]; then
    BUILD_TYPE="${2:-Release}"
    echo "Building in development mode (${BUILD_TYPE})..."
    cd "$PROJECT_ROOT"
    
    # Set build type for scikit-build-core
    export CMAKE_BUILD_TYPE="$BUILD_TYPE"
    uv pip install -e ".[dev]"
    
elif [ "$1" = "debug" ]; then
    echo "Building debug version directly with CMake..."
    BUILD_TYPE="Debug"
    BUILD_DIR="$PROJECT_ROOT/build/debug"
    
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" "$PROJECT_ROOT"
    make -j$(nproc)
    
    echo "Debug build complete. Executables are in: $BUILD_DIR/bin"
    echo "Libraries are in: $BUILD_DIR/lib"
    
# TODO: once pybindings are done, remove this section since "release" version will be 
# the default build type. Keep debug for development purposes.
elif [ "$1" = "cmake" ]; then
    BUILD_TYPE="${2:-Release}"
    echo "Building with CMake (${BUILD_TYPE})..."
    BUILD_DIR="$PROJECT_ROOT/build/cmake_$BUILD_TYPE"
    
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" "$PROJECT_ROOT"
    make -j$(nproc)
    
    echo "CMake build complete. Executables are in: $BUILD_DIR/bin"
    echo "Libraries are in: $BUILD_DIR/lib"
elif [ "$1" = "release" ]; then
    echo "Building release version directly with CMake..."
    BUILD_TYPE="Release"
    BUILD_DIR="$PROJECT_ROOT/build/release"
    
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" "$PROJECT_ROOT"
    make -j$(nproc)
    
    echo "Release build complete. Executables are in: $BUILD_DIR/bin"
    echo "Libraries are in: $BUILD_DIR/lib"
    
elif [ "$1" = "clean" ]; then
    echo "Cleaning all build directories..."
    cd "$PROJECT_ROOT"
    rm -rf build/*
    rm -rf dist/
    rm -rf *.egg-info/
    echo "All build artifacts cleaned."
    
elif [ "$1" = "clean-debug" ]; then
    echo "Cleaning debug build directory..."
    cd "$PROJECT_ROOT"
    rm -rf build/debug/
    rm -rf build/python/
    echo "Debug build artifacts cleaned."
    
elif [ "$1" = "clean-release" ]; then
    echo "Cleaning release build directory..."
    cd "$PROJECT_ROOT"
    rm -rf build/release/
    rm -rf build/python/
    echo "Release build artifacts cleaned."
    
elif [ "$1" = "clean-python" ]; then
    echo "Cleaning Python build artifacts..."
    cd "$PROJECT_ROOT"
    rm -rf build/python/
    rm -rf dist/
    rm -rf *.egg-info/
    echo "Python build artifacts cleaned."
    
else
    echo "Usage: $0 {wheel|dev|cmake|debug|release|clean|clean-debug|clean-release|clean-python} [build_type]"
    echo ""
    echo "  wheel        - Build a wheel using uv build [Release|Debug] (outputs to build/python)"
    echo "  dev          - Install in development mode using uv pip install -e . [Release|Debug] (outputs to build/python)"
    echo "  release      - CMake release build → build/release/bin/"
    echo "  debug        - CMake debug build → build/debug/bin/"
    echo "  clean        - Remove all build artifacts"
    echo "  clean-debug  - Remove only debug build artifacts"
    echo "  clean-release- Remove only release build artifacts"
    echo "  clean-python - Remove only Python build artifacts"
    echo ""
    echo "Examples:"
    echo "  $0 debug                  # Quick debug build → build/debug/bin/"
    echo "  $0 release                # Quick release build → build/release/bin/"
    echo "  $0 wheel                  # Python wheel → build/python/"
    echo "  $0 dev Debug              # Development mode → build/python/"
    echo "  $0 cmake Release          # Direct CMake → build/cmake_Release/"
    echo "  $0 clean"
    exit 1
fi
