#!/bin/bash
# Comprehensive validation script for all test cases

set -e

echo "======================================================================="
echo "  COMPREHENSIVE VALIDATION: All Test Cases"
echo "======================================================================="
echo ""

# Array of test files
TEST_FILES=(
    "trial-runList"
    "trial-input-output"
    "exclude-runList"
    "advanced-input"
    "loopCO"
)

# Results directory
RESULTS_DIR="tools/validation_results"
mkdir -p "$RESULTS_DIR"

echo "Found test files:"
for test in "${TEST_FILES[@]}"; do
    if [ -f "tools/${test}.json" ]; then
        echo "  ✓ tools/${test}.json"
    else
        echo "  ✗ tools/${test}.json (missing)"
    fi
done
echo ""

# Run C++ validation for each test
echo "Running C++ calculations..."
echo "-----------------------------------------------------------------------"

TOTAL=0
PASSED=0
FAILED=0

for test in "${TEST_FILES[@]}"; do
    JSON_FILE="tools/${test}.json"
    RESULT_FILE="$RESULTS_DIR/${test}-cpp.json"

    if [ ! -f "$JSON_FILE" ]; then
        echo "⊘ Skipping $test (JSON not found)"
        continue
    fi

    echo ""
    echo "Testing: $test"
    echo "  Input:  $JSON_FILE"
    echo "  Output: $RESULT_FILE"

    if ./tools/batch_validate "$JSON_FILE" "$RESULT_FILE" 2>&1 | grep -q "Passed:"; then
        # Extract pass/fail counts
        RESULT=$(./tools/batch_validate "$JSON_FILE" "$RESULT_FILE" 2>&1 | tail -5)
        echo "$RESULT" | grep "Passed:"
        echo "$RESULT" | grep "Failed:"

        # Check if all passed
        if echo "$RESULT" | grep -q "Failed: 0"; then
            echo "  ✓ All tests passed"
            PASSED=$((PASSED + 1))
        else
            echo "  ✗ Some tests failed"
            FAILED=$((FAILED + 1))
        fi
    else
        echo "  ✗ Validation failed"
        FAILED=$((FAILED + 1))
    fi

    TOTAL=$((TOTAL + 1))
done

echo ""
echo "======================================================================="
echo "  SUMMARY"
echo "======================================================================="
echo "Total test suites: $TOTAL"
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All validation suites passed!"
    exit 0
else
    echo "✗ $FAILED validation suite(s) failed"
    exit 1
fi
