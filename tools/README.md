# Thermochimica Validation Tools

## batch_validate

A command-line tool for running batch thermodynamic calculations and validating results against reference data.

### Usage

```bash
batch_validate <input.json> [output.json]
```

- `input.json`: JSON file containing test cases
- `output.json`: Output file for results (default: `results.json`)

### Input JSON Format

```json
{
  "test_cases": [
    {
      "name": "Test case name",
      "database": "database_file.dat",
      "temperature": 1000.0,
      "pressure": 1.0,
      "composition": {
        "Element1": mass1,
        "Element2": mass2
      },
      "reference": {
        "gibbs": -432000.0
      }
    }
  ]
}
```

**Required fields:**
- `name`: Descriptive name for the test case
- `database`: Thermodynamic database file (e.g., "CO.dat", "NobleMetals-Kaye.dat")
- `temperature`: Temperature in Kelvin
- `pressure`: Pressure in atmospheres
- `composition`: Map of element symbols to masses (grams or moles depending on units)

**Optional fields:**
- `reference`: Reference values for comparison
  - `gibbs`: Expected Gibbs energy in Joules

### Output JSON Format

```json
{
  "results": [
    {
      "name": "Test case name",
      "success": true,
      "error_code": 0,
      "gibbs_energy": -432145.67,
      "chemical_potentials": {
        "C": -123456.78,
        "O": -234567.89
      },
      "comparison": {
        "gibbs_error": 145.67
      }
    }
  ]
}
```

### Example

Run the example test suite:

```bash
cd build
./batch_validate ../tools/example_tests.json example_results.json
```

This will:
1. Read 6 test cases from `example_tests.json`
2. Run each calculation using the appropriate database and conditions
3. Write results to `example_results.json`
4. Print a summary to console

### Use Cases

#### 1. Regression Testing

Create a reference file with known-good results:
```bash
# Generate reference from current version
./batch_validate test_suite.json reference_v2.0.json

# After making changes, compare
./batch_validate test_suite.json current_results.json

# Use diff or custom script to compare JSON files
diff reference_v2.0.json current_results.json
```

#### 2. Comparing C++ vs Fortran

1. Run calculations with Fortran version, export results as JSON
2. Run same test cases with C++ version using `batch_validate`
3. Compare the output JSON files

Example comparison script:
```python
import json
import sys

def compare_results(file1, file2, tolerance=1e-6):
    with open(file1) as f1, open(file2) as f2:
        data1 = json.load(f1)
        data2 = json.load(f2)

    for r1, r2 in zip(data1['results'], data2['results']):
        name = r1['name']
        g1 = r1['gibbs_energy']
        g2 = r2['gibbs_energy']
        diff = abs(g1 - g2)

        if diff > tolerance:
            print(f"DIFF: {name}: {diff:.3e} J")
        else:
            print(f"OK: {name}")

if __name__ == "__main__":
    compare_results(sys.argv[1], sys.argv[2])
```

#### 3. Model Validation

Test each thermodynamic model with representative systems:

```json
{
  "test_cases": [
    {"name": "IDMX_test", "database": "CO.dat", ...},
    {"name": "RKMP_test", "database": "NobleMetals-Kaye.dat", ...},
    {"name": "SUBL_test", "database": "FeTiVO.dat", ...},
    {"name": "SUBG_test", "database": "CsI-Pham.dat", ...},
    {"name": "SUBQ_test", "database": "ZIRC-noSUBI.dat", ...}
  ]
}
```

#### 4. Temperature/Composition Sweeps

Generate input files programmatically to sweep parameters:

```python
import json

test_cases = []
for T in range(500, 2500, 100):
    test_cases.append({
        "name": f"MoRu_{T}K",
        "database": "NobleMetals-Kaye.dat",
        "temperature": float(T),
        "pressure": 1.0,
        "composition": {"Mo": 0.5, "Ru": 0.5}
    })

with open("temperature_sweep.json", "w") as f:
    json.dump({"test_cases": test_cases}, f, indent=2)
```

### Building

The tool is built automatically when `THERMO_BUILD_TOOLS` is ON (default):

```bash
mkdir build && cd build
cmake ..
make batch_validate
```

To disable building tools:
```bash
cmake -DTHERMO_BUILD_TOOLS=OFF ..
```

### Error Codes

Common error codes from Thermochimica:
- `0`: Success
- `10`: Data file not found
- `20`: Invalid temperature
- `21`: Invalid pressure
- `40`: No composition specified
- `-999`: Exception occurred (see error message)

See `include/thermochimica/ErrorCodes.hpp` for complete list.

### Notes

- Database files are resolved using the `DATA_DIRECTORY` CMake variable
- Standard units (K, atm, moles) are used by default
- Chemical potentials are reported in J/mol
- Gibbs energy is reported in Joules

### Future Enhancements

Potential additions:
- Support for phase amount comparison
- Statistical analysis (mean, std dev of errors)
- Parallel execution for large test suites
- CSV output format
- Plotting results (via Python script)
