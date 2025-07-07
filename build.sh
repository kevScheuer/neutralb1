#!/bin/bash

# Build script for neutralb1 project
# This script sources the GlueX environment and builds the project

set -e  # Exit on any error

echo "Setting up GlueX environment..."
source config/setup_gluex.sh

echo "Building neutralb1 project..."
make "$@"

echo "Build complete!"