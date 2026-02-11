#!/bin/bash
# Script to clean up in-source build artifacts from git tracking

set -e

echo "========================================="
echo "Cleaning up in-source build artifacts"
echo "========================================="

# Remove CMake artifacts from git tracking (but keep files on disk for now)
echo ""
echo "Removing CMake artifacts from git tracking..."

# Root directory artifacts
git rm -r --cached CMakeFiles/ 2>/dev/null || true
git rm --cached CMakeCache.txt 2>/dev/null || true
git rm --cached Makefile 2>/dev/null || true
git rm --cached cmake_install.cmake 2>/dev/null || true
git rm --cached CTestTestfile.cmake 2>/dev/null || true
git rm --cached ThermochimicaConfig.cmake 2>/dev/null || true
git rm --cached ThermochimicaConfigVersion.cmake 2>/dev/null || true

# Build executables
git rm --cached example_basic 2>/dev/null || true

# _deps directory (FetchContent downloads)
git rm -r --cached _deps/ 2>/dev/null || true

# Testing artifacts
git rm -r --cached Testing/ 2>/dev/null || true

# Legacy directory (if it's build artifacts)
git rm -r --cached legacy/ 2>/dev/null || true

# Tools directory artifacts
echo "Cleaning tools/ directory..."
git rm -r --cached tools/CMakeFiles/ 2>/dev/null || true
git rm --cached tools/Makefile 2>/dev/null || true
git rm --cached tools/cmake_install.cmake 2>/dev/null || true
git rm --cached tools/batch_validate 2>/dev/null || true
git rm --cached tools/check_elements 2>/dev/null || true
git rm --cached tools/diagnose_graphite 2>/dev/null || true
git rm --cached tools/diagnose_with_trace 2>/dev/null || true

# Tests directory artifacts
echo "Cleaning tests/ directory..."
git rm -r --cached tests/CMakeFiles/ 2>/dev/null || true
git rm --cached tests/Makefile 2>/dev/null || true
git rm --cached tests/cmake_install.cmake 2>/dev/null || true
git rm --cached tests/CTestTestfile.cmake 2>/dev/null || true
git rm --cached tests/integration_tests 2>/dev/null || true
git rm --cached tests/unit_tests 2>/dev/null || true
git rm -r --cached tests/Testing/ 2>/dev/null || true
git rm --cached tests/integration_tests[1]_include.cmake 2>/dev/null || true
git rm --cached tests/integration_tests[1]_tests.cmake 2>/dev/null || true
git rm --cached tests/unit_tests[1]_include.cmake 2>/dev/null || true
git rm --cached tests/unit_tests[1]_tests.cmake 2>/dev/null || true

echo ""
echo "✅ Build artifacts removed from git tracking"
echo ""
echo "Note: Files still exist on disk. To remove them physically:"
echo "  make clean  # or manually: rm -rf CMakeFiles/ CMakeCache.txt Makefile ..."
echo ""
echo "========================================="
echo "Next steps:"
echo "========================================="
echo "1. Review changes: git status"
echo "2. Commit the removal: git commit -m 'Remove build artifacts from git tracking'"
echo "3. Always use out-of-source builds:"
echo "     rm -rf CMakeCache.txt CMakeFiles/ Makefile"
echo "     mkdir -p build && cd build"
echo "     cmake .."
echo "     make"
echo "========================================="
