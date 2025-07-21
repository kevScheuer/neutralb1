cmake_minimum_required(VERSION 3.15)

project(neutralb1)

# This project uses Makefile for C++ builds and Python packaging via scikit-build-core
# The CMakeLists.txt is minimal since the real build system is the Makefile for now

if(SKBUILD)
    # For scikit-build-core: we're currently building a pure Python package
    # The C++ components are built separately via Makefile
    
    # For integrating C++ extensions later:
    # 1. Add custom commands to invoke your Makefile
    # 2. Use execute_process to run make commands
    # 3. Or gradually migrate Makefile rules to CMake
    
    # For now, just let scikit-build-core handle the Python package
    message(STATUS "Building Python package only - C++ components built via Makefile")
endif()
