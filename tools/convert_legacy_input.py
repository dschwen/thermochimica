#!/usr/bin/env python3
"""
Convert legacy Thermochimica input files (.ti) to JSON format for batch_validate

Usage: python convert_legacy_input.py input.ti [output.json]
"""

import sys
import json
import re
from pathlib import Path

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

def parse_legacy_input(filename):
    """Parse legacy .ti input file"""

    # Configuration
    config = {
        'temperature_unit': 'K',
        'pressure_unit': 'atm',
        'mass_unit': 'moles',
        'data_file': '',
        'nEl': 0,
        'iEl': [],
        'nCalc': 0,
        'print_mode': 0,
        'debug_mode': False
    }

    calculations = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    reading_calcs = False

    for line_num, line in enumerate(lines, 1):
        # Remove leading/trailing whitespace
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Skip comment lines (! @ # $ % & * / \ ? |)
        if line[0] in '!@#$%&*/\\?|':
            continue

        # If reading calculations (after nCalc is set)
        if reading_calcs:
            parts = line.split()
            if len(parts) >= 2 + config['nEl']:
                temp = float(parts[0])
                press = float(parts[1])
                masses = [float(parts[2 + i]) for i in range(config['nEl'])]
                calculations.append({
                    'temperature': temp,
                    'pressure': press,
                    'masses': masses
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
            # Extract just filename if it includes path
            filename = Path(value).name
            # Apply database name mapping if needed
            filename = DATABASE_MAPPING.get(filename, filename)
            config['data_file'] = filename
        elif tag in ['print mode', 'print_mode']:
            config['print_mode'] = int(value)
        elif tag in ['debug mode', 'debug_mode']:
            config['debug_mode'] = value.lower() in ['true', 't', '1', 'yes']
        elif tag in ['nel', 'nelement', 'nelements']:
            config['nEl'] = int(value)
        elif tag in ['iel', 'ielement']:
            config['iEl'] = [int(x) for x in value.split()]
        elif tag in ['ncalc', 'ncalculation', 'ncalculations']:
            config['nCalc'] = int(value)
            reading_calcs = True

    return config, calculations

def create_json_output(config, calculations, input_basename):
    """Create JSON output in batch_validate format"""

    test_cases = []

    # Get element symbols
    element_symbols = [ELEMENT_SYMBOLS[z] for z in config['iEl']]

    # Check if units are standard (K, atm, moles)
    if (config['temperature_unit'].lower() != 'k' or
        config['pressure_unit'].lower() != 'atm' or
        config['mass_unit'].lower() != 'moles'):
        print(f"Warning: Non-standard units detected:")
        print(f"  Temperature: {config['temperature_unit']} (expected K)")
        print(f"  Pressure: {config['pressure_unit']} (expected atm)")
        print(f"  Mass: {config['mass_unit']} (expected moles)")
        print(f"  Conversion may be needed!")

    for i, calc in enumerate(calculations, 1):
        # Create composition dictionary
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

    return {
        "test_cases": test_cases
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_legacy_input.py input.ti [output.json]")
        print("\nConverts legacy Thermochimica input files to JSON format")
        sys.exit(1)

    input_file = sys.argv[1]

    # Determine output filename
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # Default: replace .ti with .json
        output_file = Path(input_file).stem + '.json'

    print(f"Converting {input_file} -> {output_file}")

    # Parse input
    config, calculations = parse_legacy_input(input_file)

    print(f"\nParsed configuration:")
    print(f"  Database: {config['data_file']}")
    print(f"  Elements: {[ELEMENT_SYMBOLS[z] for z in config['iEl']]} (Z = {config['iEl']})")
    print(f"  Units: {config['temperature_unit']}, {config['pressure_unit']}, {config['mass_unit']}")
    print(f"  Calculations: {len(calculations)}")

    # Create JSON
    input_basename = Path(input_file).stem
    json_data = create_json_output(config, calculations, input_basename)

    # Write output
    with open(output_file, 'w') as f:
        json.dump(json_data, f, indent=2)

    print(f"\nConversion complete!")
    print(f"Created {len(json_data['test_cases'])} test cases")
    print(f"\nRun with: batch_validate {output_file} results.json")

if __name__ == '__main__':
    main()
