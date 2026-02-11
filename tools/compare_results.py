#!/usr/bin/env python3
"""
Compare C++ and Fortran Thermochimica results

Usage: python compare_results.py cpp_results.json fortran_results.json [tolerance]
"""

import json
import sys

def compare_results(cpp_file, fortran_file, tolerance=1e-3):
    """Compare C++ and Fortran results"""

    with open(cpp_file) as f:
        cpp_data = json.load(f)

    with open(fortran_file) as f:
        fortran_data = json.load(f)

    cpp_results = cpp_data['results']
    fortran_results = fortran_data.get('results', fortran_data)

    # Handle case where fortran output uses numbered keys ("1", "2", "3", ...)
    if isinstance(fortran_results, dict) and '1' in fortran_results:
        # Convert numbered dict to list
        fortran_list = []
        for i in range(1, len(cpp_results) + 1):
            if str(i) in fortran_results:
                fortran_list.append(fortran_results[str(i)])
        fortran_results = fortran_list
    elif isinstance(fortran_results, dict) and 'results' in fortran_results:
        fortran_results = fortran_results['results']

    print(f"Comparing {len(cpp_results)} test cases...")
    print(f"Tolerance: {tolerance} J")
    print(f"{'='*70}\n")

    passed = 0
    failed = 0
    max_diff = 0.0
    max_diff_case = ""

    for i, cpp in enumerate(cpp_results):
        name = cpp['name']

        # Try to match by name or index
        if isinstance(fortran_results, list) and i < len(fortran_results):
            fortran = fortran_results[i]
        elif isinstance(fortran_results, dict):
            fortran = fortran_results.get(name, {})
        else:
            fortran = {}

        if not cpp['success']:
            print(f"❌ {name}: C++ calculation failed (code {cpp['error_code']})")
            failed += 1
            continue

        if not fortran.get('success', True):
            print(f"❌ {name}: Fortran calculation failed")
            failed += 1
            continue

        cpp_g = cpp['gibbs_energy']
        # Handle different field names in Fortran output
        fortran_g = fortran.get('gibbs_energy',
                                 fortran.get('gibbs',
                                 fortran.get('integral Gibbs energy', 0.0)))

        diff = abs(cpp_g - fortran_g)
        rel_diff = abs(diff / fortran_g) if abs(fortran_g) > 1e-10 else 0

        if diff > max_diff:
            max_diff = diff
            max_diff_case = name

        if diff > tolerance:
            print(f"❌ {name}:")
            print(f"   C++:     {cpp_g:>16.6e} J")
            print(f"   Fortran: {fortran_g:>16.6e} J")
            print(f"   Diff:    {diff:>16.6e} J ({rel_diff*100:>6.3f}%)\n")
            failed += 1
        else:
            print(f"✓ {name:30s} Diff: {diff:>12.3e} J")
            passed += 1

    print(f"\n{'='*70}")
    print(f"Results Summary:")
    print(f"  Passed: {passed}/{len(cpp_results)}")
    print(f"  Failed: {failed}/{len(cpp_results)}")
    print(f"  Max difference: {max_diff:.3e} J (in {max_diff_case})")

    if passed == len(cpp_results):
        print(f"\n✓ All tests passed!")
        return True
    else:
        print(f"\n❌ {failed} tests failed")
        return False

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_results.py cpp_results.json fortran_results.json [tolerance]")
        print("\nCompares C++ and Fortran Thermochimica results")
        print("Default tolerance: 1e-3 J")
        sys.exit(1)

    cpp_file = sys.argv[1]
    fortran_file = sys.argv[2]
    tolerance = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3

    try:
        success = compare_results(cpp_file, fortran_file, tolerance)
        sys.exit(0 if success else 1)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(2)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON - {e}")
        sys.exit(2)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(2)

if __name__ == '__main__':
    main()
