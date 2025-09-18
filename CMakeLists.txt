cmake_minimum_required(VERSION 3.15)

project(neutralb1 VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable position independent code for Python bindings
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Project directories - set early so they can be used throughout
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
set(SRC_DIR "${PROJECT_ROOT}/src/neutralb1_cpp")
set(INCLUDE_DIR "${PROJECT_ROOT}/include")

# Build type handling
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Normalize build type case (debug or release)
string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE_LOWER)
if(BUILD_TYPE_LOWER STREQUAL "debug")
    set(CMAKE_BUILD_TYPE Debug)
elseif(BUILD_TYPE_LOWER STREQUAL "release") 
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags based on build type
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -DDEBUG")
    message(STATUS "Building Debug version with -g -O0")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG") 
    message(STATUS "Building Release version with -O2")
endif()

# Set output directories for direct cmake builds
if(NOT SKBUILD)
    string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE_LOWER)
    
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_ROOT}/build/${BUILD_TYPE_LOWER}/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_ROOT}/build/${BUILD_TYPE_LOWER}/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_ROOT}/build/${BUILD_TYPE_LOWER}/lib")
    
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_ROOT}/build/debug/bin")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_ROOT}/build/release/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG "${PROJECT_ROOT}/build/debug/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE "${PROJECT_ROOT}/build/release/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${PROJECT_ROOT}/build/debug/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${PROJECT_ROOT}/build/release/lib")

    # Create the directories
    file(MAKE_DIRECTORY "${PROJECT_ROOT}/build/${BUILD_TYPE_LOWER}/bin")
    file(MAKE_DIRECTORY "${PROJECT_ROOT}/build/${BUILD_TYPE_LOWER}/lib")
    
    message(STATUS "Output directories set to: ${PROJECT_ROOT}/build/${BUILD_TYPE_LOWER}/")
endif()

# Find required packages
find_package(PkgConfig REQUIRED)

# Environment detection - try to use environment variables from setup script
if(DEFINED ENV{ROOTSYS})
    message(STATUS "Using ROOT from environment: $ENV{ROOTSYS}")
    set(ROOT_DIR "$ENV{ROOTSYS}")
    
    # Set ROOT variables manually if not found by find_package
    set(ROOT_INCLUDE_DIRS "$ENV{ROOTSYS}/include")
    set(ROOT_LIBRARY_DIR "$ENV{ROOTSYS}/lib")
    
    # Get ROOT libraries
    execute_process(
        COMMAND $ENV{ROOTSYS}/bin/root-config --libs
        OUTPUT_VARIABLE ROOT_LIBRARIES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    # Get ROOT compile flags
    execute_process(
        COMMAND $ENV{ROOTSYS}/bin/root-config --cflags
        OUTPUT_VARIABLE ROOT_CFLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    # Parse ROOT flags and add them to compile options
    string(REPLACE " " ";" ROOT_CFLAGS_LIST ${ROOT_CFLAGS})
    
else()
    # Fallback: try to find ROOT normally
    find_package(ROOT REQUIRED)
    if(ROOT_FOUND)
        set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS})
        set(ROOT_LIBRARIES ${ROOT_LIBRARIES})
    else()
        message(FATAL_ERROR "ROOT not found. Please source your setup script first.")
    endif()
endif()

# Include the main include directory and all subdirectories
include_directories("${INCLUDE_DIR}")
file(GLOB_RECURSE INCLUDE_SUBDIRS "${INCLUDE_DIR}/*")
foreach(subdir ${INCLUDE_SUBDIRS})
    if(IS_DIRECTORY "${subdir}")
        include_directories("${subdir}")
    endif()
endforeach()
include_directories(${ROOT_INCLUDE_DIRS})

# Add environment-based include directories if available
if(DEFINED ENV{AMPTOOLS})
    include_directories("$ENV{AMPTOOLS}")
    link_directories("$ENV{AMPTOOLS}/lib")
endif()

if(DEFINED ENV{HALLD_SIM_HOME})
    include_directories("$ENV{HALLD_SIM_HOME}/src/libraries")
    include_directories("$ENV{HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/include")
    link_directories("$ENV{HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib")
endif()

# Add ROOT library directory
if(ROOT_LIBRARY_DIR)
    link_directories(${ROOT_LIBRARY_DIR})
endif()

# Collect source files
file(GLOB_RECURSE BATCH_SOURCES "${SRC_DIR}/batch/*.cc")
file(GLOB_RECURSE ANALYSIS_SOURCES "${SRC_DIR}/analysis/*.cc")
file(GLOB_RECURSE SELECTION_SOURCES "${SRC_DIR}/selection/*.cc")
file(GLOB_RECURSE UTILS_SOURCES "${SRC_DIR}/utils/*.cc")

# Remove shell scripts and ROOT macros from executable list
list(FILTER BATCH_SOURCES EXCLUDE REGEX "\\.(sh|C)$")

# Create a static library from utility, analysis, and selection code
if(UTILS_SOURCES OR ANALYSIS_SOURCES OR SELECTION_SOURCES)
    add_library(neutralb1_core STATIC
        ${UTILS_SOURCES} 
        ${ANALYSIS_SOURCES} 
        ${SELECTION_SOURCES}
    )
    
    # Link ROOT libraries to the core library
    target_link_libraries(neutralb1_core 
        ${ROOT_LIBRARIES}
    )
    
    # Add environment-specific libraries if available
    if(DEFINED ENV{AMPTOOLS})
        target_link_libraries(neutralb1_core AmpTools AMPTOOLS_AMPS)
    endif()
    
    # Apply ROOT compile flags
    if(ROOT_CFLAGS_LIST)
        target_compile_options(neutralb1_core PRIVATE ${ROOT_CFLAGS_LIST})
    endif()
endif()

# Create executables from batch sources
foreach(BATCH_SOURCE ${BATCH_SOURCES})
    # Get filename without extension
    get_filename_component(BATCH_NAME ${BATCH_SOURCE} NAME_WE)
    
    # Create executable
    add_executable(${BATCH_NAME} ${BATCH_SOURCE})
    
    # Link with core library and ROOT
    if(TARGET neutralb1_core)
        target_link_libraries(${BATCH_NAME} neutralb1_core)
    endif()
    
    target_link_libraries(${BATCH_NAME} 
        ${ROOT_LIBRARIES}
    )
    
    # Add environment-specific libraries if available
    if(DEFINED ENV{AMPTOOLS})
        target_link_libraries(${BATCH_NAME} AmpTools AMPTOOLS_AMPS)
    endif()
    
    # Apply ROOT compile flags
    if(ROOT_CFLAGS_LIST)
        target_compile_options(${BATCH_NAME} PRIVATE ${ROOT_CFLAGS_LIST})
    endif()
endforeach()

# Installation rules
if(SKBUILD)
    # When building with scikit-build-core, install to wheel
    message(STATUS "Building with scikit-build-core for ${CMAKE_BUILD_TYPE}")
    
    # Install the core library
    if(TARGET neutralb1_core)
        install(TARGETS neutralb1_core
            LIBRARY DESTINATION neutralb1
            RUNTIME DESTINATION neutralb1
        )
    endif()
    
    # TODO: replace this with pybindings
    # foreach(BATCH_SOURCE ${BATCH_SOURCES})
    #     get_filename_component(BATCH_NAME ${BATCH_SOURCE} NAME_WE)
    #     if(TARGET ${BATCH_NAME})
    #         install(TARGETS ${BATCH_NAME}
    #             RUNTIME DESTINATION neutralb1/bin
    #         )
    #     endif()
    # endforeach()
    
else()
    # Normal CMake build
    string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE_LOWER)
    
    # Installation rules (optional, since binaries are already in correct location)
    set(CMAKE_INSTALL_PREFIX "${PROJECT_ROOT}/build/${BUILD_TYPE_LOWER}" CACHE PATH "Install prefix")
    
    if(TARGET neutralb1_core)
        install(TARGETS neutralb1_core
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
        )
    endif()
    
    foreach(BATCH_SOURCE ${BATCH_SOURCES})
        get_filename_component(BATCH_NAME ${BATCH_SOURCE} NAME_WE)
        if(TARGET ${BATCH_NAME})
            install(TARGETS ${BATCH_NAME}
                RUNTIME DESTINATION bin
            )
        endif()
    endforeach()
    
    # Install headers
    install(DIRECTORY ${INCLUDE_DIR}/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h"
    )
endif()

# Custom targets
add_custom_target(info
    COMMAND ${CMAKE_COMMAND} -E echo "Project Root: ${PROJECT_ROOT}"
    COMMAND ${CMAKE_COMMAND} -E echo "Build Type: ${CMAKE_BUILD_TYPE}"
    COMMAND ${CMAKE_COMMAND} -E echo "Install Prefix: ${CMAKE_INSTALL_PREFIX}"
    COMMAND ${CMAKE_COMMAND} -E echo "ROOT Directory: ${ROOT_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "Source files found:"
    COMMAND ${CMAKE_COMMAND} -E echo "  Batch: ${BATCH_SOURCES}"
    COMMAND ${CMAKE_COMMAND} -E echo "  Analysis: ${ANALYSIS_SOURCES}"
    COMMAND ${CMAKE_COMMAND} -E echo "  Selection: ${SELECTION_SOURCES}"
    COMMAND ${CMAKE_COMMAND} -E echo "  Utils: ${UTILS_SOURCES}"
    COMMAND ${CMAKE_COMMAND} -E echo "Compiler: ${CMAKE_CXX_COMPILER}"
    COMMAND ${CMAKE_COMMAND} -E echo "Flags: ${CMAKE_CXX_FLAGS}"
)

# Clean targets matching Makefile functionality
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning all build directories..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/build"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_ROOT}/build"
    COMMENT "Clean all build artifacts"
)

add_custom_target(clean-debug
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning debug build directory..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/build/debug"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/build/python"
    COMMENT "Clean debug build artifacts"
)

add_custom_target(clean-release
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning release build directory..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/build/release"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/build/python"
    COMMENT "Clean release build artifacts"
)

add_custom_target(clean-python
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning Python build artifacts..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/build/python"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/dist"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PROJECT_ROOT}/*.egg-info"
    COMMENT "Clean Python build artifacts"
)

# Environment setup helper
add_custom_target(check-env
    COMMAND ${CMAKE_COMMAND} -E echo "Checking environment variables..."
    COMMAND ${CMAKE_COMMAND} -E echo "ROOTSYS: $ENV{ROOTSYS}"
    COMMAND ${CMAKE_COMMAND} -E echo "AMPTOOLS: $ENV{AMPTOOLS}"
    COMMAND ${CMAKE_COMMAND} -E echo "HALLD_SIM_HOME: $ENV{HALLD_SIM_HOME}"
)
