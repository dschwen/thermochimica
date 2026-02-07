#!/usr/bin/env python3
"""
Convert legacy Thermochimica input files (.ti) to JSON format for batch_validate

Supports both explicit calculations and range/loop notation:
- Explicit: nCalc=3 followed by explicit T P mass1 mass2 lines
- Range: temperature = 900:1000:50 (generates 3 points: 900, 950, 1000)
- Range: mass(6) = 1:3:1 (generates 3 points: 1, 2, 3)

Usage: python convert_legacy_input_v2.py input.ti [output.json]
"""

import sys
import json
import re
import numpy as np
from pathlib import Path
from itertools import product

# Atomic number to element symbol mapping
ELEMENT_SYMBOLS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
    19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe',
    27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se',
    35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
    43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce',
    59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy',
    67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W',
    75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
    83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf',
    99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr'
}

# Database filename mapping (legacy -> current)
DATABASE_MAPPING = {
    'C-O.dat': 'CO.dat',
    'Kaye_NobleMetals.dat': 'NobleMetals-Kaye.dat',
}

def parse_range(value_str):
    """
    Parse range notation: start:end[:step]

    Returns tuple (is_range, values)
    - is_range: bool indicating if this is a range
    - values: list of values or single value

    Examples:
        "1.0" -> (False, 1.0)
        "900:1000:50" -> (True, [900, 950, 1000])
        "1:5" -> (True, [1, 2, 3, 4, 5])
        "10:5:-1" -> (True, [10, 9, 8, 7, 6, 5])
    """
    if ':' not in value_str:
        # Not a range, just a single value
        try:
            return (False, float(value_str))
        except ValueError:
            return (False, value_str)

    # Parse range notation
    parts = value_str.split(':')

    if len(parts) == 2:
        # start:end (default step = 1)
        start = float(parts[0])
        end = float(parts[1])
        step = 1.0 if end >= start else -1.0
    elif len(parts) == 3:
        # start:end:step
        start = float(parts[0])
        end = float(parts[1])
        step = float(parts[2])
    else:
        raise ValueError(f"Invalid range format: {value_str}")

    # Generate values using numpy
    # Use arange but ensure endpoint is included if it's exactly on a step
    if step > 0:
        values = np.arange(start, end + step/2, step)
    else:
        values = np.arange(start, end + step/2, step)

    return (True, list(values))

def parse_legacy_input(filename):
    """Parse legacy .ti input file (both explicit and range formats)"""

    # Configuration
    config = {
        'temperature_unit': 'K',
        'pressure_unit': 'atm',
        'mass_unit': 'moles',
        'data_file': '',
        'print_mode': 0,
        'debug_mode': False,
        'step_together': False
    }

    # Variables that can be ranges or scalars
    temperature = None
    pressure = None
    masses = {}  # atomic_number -> value(s)

    # For explicit calculation format
    explicit_calculations = []
    elements_list = []  # iEl list
    reading_explicit_calcs = False
    nCalc = 0

    with open(filename, 'r') as f:
        lines = f.readlines()

    for line_num, line in enumerate(lines, 1):
        line = line.strip()

        # Skip empty lines and comments
        if not line or line[0] in '!@#$%&*/\\?|':
            continue

        # If reading explicit calculations
        if reading_explicit_calcs and len(explicit_calculations) < nCalc:
            parts = line.split()
            if len(parts) >= 2 + len(elements_list):
                temp = float(parts[0])
                press = float(parts[1])
                calc_masses = [float(parts[2 + i]) for i in range(len(elements_list))]
                explicit_calculations.append({
                    'temperature': temp,
                    'pressure': press,
                    'masses': calc_masses
                })
            continue

        # Look for = delimiter
        if '=' not in line:
            continue

        # Split on =
        parts = line.split('=', 1)
        tag = parts[0].strip().lower()
        value = parts[1].strip()

        # Parse based on tag
        if tag in ['temperature unit', 't_unit', 'temperature_unit']:
            config['temperature_unit'] = value
        elif tag in ['pressure unit', 'p_unit', 'pressure_unit']:
            config['pressure_unit'] = value
        elif tag in ['mass unit', 'm_unit', 'mass_unit']:
            config['mass_unit'] = value
        elif tag in ['data file', 'data', 'data_file', 'dat']:
            filename = Path(value).name
            filename = DATABASE_MAPPING.get(filename, filename)
            config['data_file'] = filename
        elif tag in ['print mode', 'print_mode']:
            config['print_mode'] = int(value)
        elif tag in ['debug mode', 'debug_mode']:
            config['debug_mode'] = value.lower() in ['true', 't', '1', 'yes']
        elif tag in ['step together', 'step_together']:
            config['step_together'] = value.lower() in ['true', 't', '1', 'yes']
        elif tag in ['temperature', 'temp', 't']:
            temperature = parse_range(value)
        elif tag in ['pressure', 'press', 'p']:
            pressure = parse_range(value)
        elif tag in ['nel', 'nelement', 'nelements']:
            # Not needed for range format
            pass
        elif tag in ['iel', 'ielement']:
            elements_list = [int(x) for x in value.split()]
        elif tag in ['ncalc', 'ncalculation', 'ncalculations']:
            nCalc = int(value)
            reading_explicit_calcs = True
        elif tag.startswith('mass(') or tag.startswith('m('):
            # Extract atomic number: mass(6) or m(6)
            match = re.search(r'\((\d+)\)', tag)
            if match:
                atomic_num = int(match.group(1))
                masses[atomic_num] = parse_range(value)

    # Return appropriate format
    if explicit_calculations:
        # Explicit format
        return 'explicit', config, elements_list, explicit_calculations
    else:
        # Range format
        return 'range', config, temperature, pressure, masses

def generate_test_cases_from_ranges(config, temperature, pressure, masses, input_basename):
    """Generate test cases from range notation"""

    # Check if we have temperature
    if temperature is None:
        print("Warning: No temperature specified!")
        return []

    # Check if we have pressure
    if pressure is None:
        print("Warning: No pressure specified, using default 1.0 atm")
        pressure = (False, 1.0)

    # Build parameter lists
    is_temp_range, temp_values = temperature
    is_press_range, press_values = pressure

    # Convert scalar to list
    if not is_temp_range:
        temp_values = [temp_values]
    if not is_press_range:
        press_values = [press_values]

    # Build mass parameter lists
    mass_params = {}
    atomic_numbers = sorted(masses.keys())

    for z in atomic_numbers:
        is_range, values = masses[z]
        if not is_range:
            values = [values]
        mass_params[z] = values

    # Generate all combinations
    # Create parameter grid
    param_lists = [temp_values, press_values]
    param_names = ['temperature', 'pressure']

    for z in atomic_numbers:
        param_lists.append(mass_params[z])
        param_names.append(f'mass_{z}')

    # Generate Cartesian product
    test_cases = []

    for i, combination in enumerate(product(*param_lists), 1):
        # Extract parameters
        temp = combination[0]
        press = combination[1]
        calc_masses = combination[2:]

        # Create composition dictionary
        composition = {}
        for j, z in enumerate(atomic_numbers):
            if z in ELEMENT_SYMBOLS:
                composition[ELEMENT_SYMBOLS[z]] = calc_masses[j]

        test_case = {
            "name": f"{input_basename}_{i:03d}",
            "database": config['data_file'],
            "temperature": temp,
            "pressure": press,
            "composition": composition
        }

        test_cases.append(test_case)

    return test_cases

def generate_test_cases_from_explicit(config, elements_list, calculations, input_basename):
    """Generate test cases from explicit calculation format"""

    test_cases = []
    element_symbols = [ELEMENT_SYMBOLS[z] for z in elements_list]

    for i, calc in enumerate(calculations, 1):
        composition = {}
        for j, symbol in enumerate(element_symbols):
            composition[symbol] = calc['masses'][j]

        test_case = {
            "name": f"{input_basename}_{i:03d}",
            "database": config['data_file'],
            "temperature": calc['temperature'],
            "pressure": calc['pressure'],
            "composition": composition
        }

        test_cases.append(test_case)

    return test_cases

def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_legacy_input_v2.py input.ti [output.json]")
        print("\nConverts legacy Thermochimica input files to JSON format")
        print("Supports both explicit calculations and range/loop notation")
        sys.exit(1)

    input_file = sys.argv[1]

    # Determine output filename
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        output_file = Path(input_file).stem + '.json'

    print(f"Converting {input_file} -> {output_file}")

    # Parse input
    result = parse_legacy_input(input_file)
    format_type = result[0]

    if format_type == 'explicit':
        _, config, elements_list, calculations = result
        print(f"\nFormat: Explicit calculations")
        print(f"  Database: {config['data_file']}")
        print(f"  Elements: {[ELEMENT_SYMBOLS[z] for z in elements_list]} (Z = {elements_list})")
        print(f"  Calculations: {len(calculations)}")

        input_basename = Path(input_file).stem
        test_cases = generate_test_cases_from_explicit(config, elements_list, calculations, input_basename)

    elif format_type == 'range':
        _, config, temperature, pressure, masses = result
        print(f"\nFormat: Range/Loop notation")
        print(f"  Database: {config['data_file']}")

        # Show ranges
        if temperature:
            is_range, values = temperature
            if is_range:
                print(f"  Temperature: {values[0]:.1f} to {values[-1]:.1f} ({len(values)} points)")
            else:
                print(f"  Temperature: {values} K")

        if pressure:
            is_range, values = pressure
            if is_range:
                print(f"  Pressure: {values[0]:.1f} to {values[-1]:.1f} ({len(values)} points)")
            else:
                print(f"  Pressure: {values} atm")

        print(f"  Elements with mass specified:")
        for z in sorted(masses.keys()):
            is_range, values = masses[z]
            symbol = ELEMENT_SYMBOLS.get(z, f"Z={z}")
            if is_range:
                print(f"    {symbol}: {values[0]:.2f} to {values[-1]:.2f} ({len(values)} points)")
            else:
                print(f"    {symbol}: {values:.2f}")

        input_basename = Path(input_file).stem
        test_cases = generate_test_cases_from_ranges(config, temperature, pressure, masses, input_basename)

        total_combinations = len(test_cases)
        print(f"  Total combinations: {total_combinations}")

    # Check units
    if (config['temperature_unit'].lower() != 'k' or
        config['pressure_unit'].lower() != 'atm' or
        config['mass_unit'].lower() != 'moles'):
        print(f"\nWarning: Non-standard units detected:")
        print(f"  Temperature: {config['temperature_unit']} (expected K)")
        print(f"  Pressure: {config['pressure_unit']} (expected atm)")
        print(f"  Mass: {config['mass_unit']} (expected moles)")

    # Create JSON output
    json_data = {"test_cases": test_cases}

    # Write output
    with open(output_file, 'w') as f:
        json.dump(json_data, f, indent=2)

    print(f"\nConversion complete!")
    print(f"Created {len(test_cases)} test cases")
    print(f"\nRun with: ./tools/batch_validate {output_file} results.json")

if __name__ == '__main__':
    main()
