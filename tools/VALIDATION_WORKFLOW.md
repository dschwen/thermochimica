# Validation Workflow: Comparing C++ with Fortran

This guide explains how to validate the C++ implementation against the legacy Fortran version.

## Overview

The validation workflow consists of:
1. Convert legacy Fortran input files (`.ti`) to JSON format
2. Run calculations using C++ `batch_validate`
3. Run calculations using Fortran `RunCalculationList`
4. Compare the JSON outputs

## Tools

### 1. convert_legacy_input.py

Converts legacy `.ti` format to JSON format compatible with `batch_validate`.

**Usage:**
```bash
python3 tools/convert_legacy_input.py inputs/trial-runList.ti tools/trial-runList.json
```

**Features:**
- Parses legacy key=value format
- Maps atomic numbers to element symbols
- Handles database filename changes (e.g., `C-O.dat` → `CO.dat`)
- Creates test cases with unique names

**Input format (.ti):**
```
temperature unit  = K
pressure unit     = atm
mass unit         = moles
data file         = data/CO.dat
nEl         = 2
iEl         = 6 8
nCalc       = 3
900  1 1 1
950  2 1 2
1000 3 2 1
```

**Output format (JSON):**
```json
{
  "test_cases": [
    {
      "name": "trial-runList_001",
      "database": "CO.dat",
      "temperature": 900.0,
      "pressure": 1.0,
      "composition": {
        "C": 1.0,
        "O": 1.0
      }
    }
  ]
}
```

### 2. batch_validate

Runs Thermochimica calculations for multiple test cases and outputs results in JSON.

**Usage:**
```bash
./tools/batch_validate tools/trial-runList.json tools/cpp-results.json
```

**Output format:**
```json
{
  "results": [
    {
      "name": "trial-runList_001",
      "success": true,
      "error_code": 0,
      "gibbs_energy": -450244.0,
      "chemical_potentials": {
        "C": 9.299,
        "O": -17.056
      }
    }
  ]
}
```

## Validation Workflow

### Step 1: Convert Legacy Inputs

Convert all legacy input files to JSON:

```bash
# Convert a single file
python3 tools/convert_legacy_input.py inputs/trial-runList.ti tools/trial-runList.json

# Convert all files (batch)
for file in inputs/*.ti; do
    basename=$(basename "$file" .ti)
    python3 tools/convert_legacy_input.py "$file" "tools/${basename}.json"
done
```

### Step 2: Run C++ Calculations

Run batch validation using the C++ implementation:

```bash
./tools/batch_validate tools/trial-runList.json tools/cpp-results.json
```

### Step 3: Run Fortran Calculations

Run the same test cases using the Fortran version:

```bash
cd legacy
./bin/RunCalculationList < ../inputs/trial-runList.ti
# Results are written to outputs/thermochimica_out.json
```

### Step 4: Compare Results

Compare the C++ and Fortran outputs using the comparison script:

```bash
python3 tools/compare_results.py tools/cpp-results.json outputs/thermochimica_out.json
```

**Manual comparison (Python):**

```python
import json
import sys

def compare_results(cpp_file, fortran_file, tolerance=1e-3):
    """Compare C++ and Fortran results"""

    with open(cpp_file) as f:
        cpp_data = json.load(f)

    with open(fortran_file) as f:
        fortran_data = json.load(f)

    cpp_results = cpp_data['results']
    fortran_results = fortran_data['results']

    print(f"Comparing {len(cpp_results)} test cases...")
    print(f"Tolerance: {tolerance} J\n")

    passed = 0
    failed = 0

    for i, (cpp, fortran) in enumerate(zip(cpp_results, fortran_results)):
        name = cpp['name']

        if not cpp['success']:
            print(f"❌ {name}: C++ calculation failed")
            failed += 1
            continue

        if not fortran.get('success', True):
            print(f"❌ {name}: Fortran calculation failed")
            failed += 1
            continue

        cpp_g = cpp['gibbs_energy']
        fortran_g = fortran.get('gibbs_energy', 0.0)

        diff = abs(cpp_g - fortran_g)
        rel_diff = abs(diff / fortran_g) if fortran_g != 0 else 0

        if diff > tolerance:
            print(f"❌ {name}:")
            print(f"   C++:     {cpp_g:.6e} J")
            print(f"   Fortran: {fortran_g:.6e} J")
            print(f"   Diff:    {diff:.6e} J ({rel_diff*100:.3f}%)")
            failed += 1
        else:
            print(f"✓ {name}: {diff:.3e} J")
            passed += 1

    print(f"\n{'='*50}")
    print(f"Passed: {passed}/{len(cpp_results)}")
    print(f"Failed: {failed}/{len(cpp_results)}")

    return failed == 0

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python compare_results.py cpp_results.json fortran_results.json")
        sys.exit(1)

    success = compare_results(sys.argv[1], sys.argv[2])
    sys.exit(0 if success else 1)
```

Save as `tools/compare_results.py` and run:

```bash
python3 tools/compare_results.py tools/cpp-results.json tools/fortran-results.json
```

## Converted Test Cases

Currently converted legacy inputs:

| File | Test Cases | Status |
|------|-----------|--------|
| `trial-runList.json` | 3 | ✓ Validated |
| `trial-input-output.json` | 3 | ✓ Validated |
| `exclude-runList.json` | 1 | ✓ Validated |
| `example_tests.json` | 6 | ✓ Validated (manual) |

## Database Name Mapping

Some legacy database names have changed:

| Legacy Name | Current Name |
|------------|--------------|
| `C-O.dat` | `CO.dat` |
| `Kaye_NobleMetals.dat` | `NobleMetals-Kaye.dat` |

These mappings are handled automatically by `convert_legacy_input.py`.

## Troubleshooting

### "Data file not found" error

- Check that the database file exists in `data/` directory
- Verify the mapping in `convert_legacy_input.py`
- Add new mappings if needed

### Parser errors

- Ensure the `.ti` file follows the expected format
- Check for special characters in values
- Verify element atomic numbers are valid (1-103)

### Result mismatches

- Check units (K, atm, moles are standard)
- Verify element composition
- Compare intermediate values (chemical potentials, phase amounts)
- Check for numerical precision differences

## Future Enhancements

- Automated comparison with statistical analysis
- Phase amount comparison
- Chemical potential comparison
- Parallel execution for large test suites
- CI/CD integration for regression testing
