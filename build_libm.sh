#!/bin/bash
# Script to build glibc's libm with AVX/FMA optimizations.
# This script configures and builds glibc in a separate directory.

set -euo pipefail

# Configuration
GLIBC_SRC="/tmp/maths/glibc"
BUILD_DIR="/tmp/maths/glibc-build"
INSTALL_DIR="/tmp/maths/glibc-install"
LOG_FILE="/tmp/maths/glibc-build.log"

# Check for source directory
if [[ ! -d "$GLIBC_SRC" ]]; then
    echo "Error: glibc source directory not found at $GLIBC_SRC"
    exit 1
fi

# Create build and install directories
mkdir -p "$BUILD_DIR"
mkdir -p "$INSTALL_DIR"
mkdir -p "$(dirname "$LOG_FILE")"

echo "Building glibc for libm optimizations..."
echo "Source: $GLIBC_SRC"
echo "Build:  $BUILD_DIR"
echo "Log:    $LOG_FILE"

cd "$BUILD_DIR"

# Configure glibc
# --enable-multiarch: Enable multiple architecture-specific implementations (AVX, FMA, etc.)
# CFLAGS="-O3 -march=native": Request maximum optimization for the current CPU
# --disable-werror: Avoid stopping on compiler warnings (common when using -O3 or newer compilers)
if [[ ! -f "config.status" ]]; then
    echo "Configuring glibc..."
    "$GLIBC_SRC/configure" \
        --prefix="$INSTALL_DIR" \
        --enable-multiarch \
        --disable-profile \
        --disable-werror \
        CFLAGS="-O3 -march=native" \
        ASFLAGS="-march=native" \
        &> "$LOG_FILE"
fi

echo "Compiling glibc (this may take a while)..."
# We build the entire glibc to ensure libm.so is correctly linked.
make -j$(nproc) &>> "$LOG_FILE"

echo "Build finished successfully."
echo "Optimized libm.so can be found at: $BUILD_DIR/math/libm.so"

# Verify built libm.so
if [[ -f "$BUILD_DIR/math/libm.so" ]]; then
    echo "Verifying symbols in built libm.so:"
    nm -D "$BUILD_DIR/math/libm.so" | grep -E " (exp|log|sin|cos)$" || true
fi
