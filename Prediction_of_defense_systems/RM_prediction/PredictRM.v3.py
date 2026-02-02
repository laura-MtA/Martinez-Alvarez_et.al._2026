#!/usr/bin/env python3
"""
PredictRM.v3.py

Purpose
-------
Classify candidate restriction-modification (RM) systems based on REBASE/MMseqs2 component hits and local genomic context (±5 proteins).

This script is used after generating:
- ar53.prot2RM.nodef.txt    (all proteins + RM annotations; sorted)
- neighborhood.txt          (subset with ±5 context around RM candidate proteins)

Input
-----
A single argument: neighborhood.txt

neighborhood.txt is expected to be a tab-delimited file with a header line and at least the columns:
- genome.acc
- protein.id
- module
- type

Other columns may be present (e.g. ref, eval, bitscore) but are not used for system calling.

Interpretation rules (as implemented here)
-----------------------------------------
- Type IV and Type IIG:
  If a row has type == "IV" or type == "IIG", it is reported directly as a predicted system.

- Type I, Type II and Type III:
  The script scans rows where type is one of: IR, IIR, IIIR
  For each such row, it collects other rows in a ±5 window around it (by row index in the input file), then checks for the presence of the required component types:

  IR   requires: IR + IS + IM
  IIR  requires: IIR + IIM
  IIIR requires: IIIR + IIIM

Output
------
results.tsv (written in the current directory)

The output is a plain text file structured as blocks separated by dashed lines, with each component line as:
  <type> <tab> <module> <tab> <genome.acc> <tab> <protein.id>

Notes
-----
- This script assumes that neighborhood.txt already represents local genomic context, produced by extracting ±5 rows around candidate proteins from an ordered protein list.
- The script prints the header and the first parsed row to stdout (debugging output).

"""

import sys

def parse_file(filename):
    """
    Read neighborhood.txt into a list of dicts (one dict per line),
    using the first line as a tab-delimited header.
    """
    data = []
    with open(filename, "r") as file:
        headers = file.readline().strip().split("\t")
        print(f"Headers: {headers}")  # Debugging step
        for line in file:
            values = line.strip().split("\t")
            row = dict(zip(headers, values))
            data.append(row)
            if len(data) == 1:  # Print the first row for debugging
                print(f"First row: {row}")  # Debugging step
    return data

def find_defense_systems(data):
       results = []
    for index, row in enumerate(data):
        if "type" in row:
            system_type = row["type"]
            genome = row["genome.acc"]
            protein = row["protein.id"]

            # Type IV and Type IIG are reported directly (no component combination required)
            if system_type in ["IV", "IIG"]:
                results.append([row])

            # Type I, II and III are evaluated based on presence of required components
            elif system_type in ["IR", "IIR", "IIIR"]:
                # Collect rows in a ±5 window around the current row (by row index in the file)
                start_index = max(0, index - 5)
                end_index = min(len(data), index + 6)

                components = [row]
                for i in range(start_index, end_index):
                    if i != index and "type" in data[i]:
                        components.append(data[i])

                # Keep only cases where the required component types are present
                if check_components(system_type, components):
                    results.append(components)

    return results


def check_components(system_type, components):
    types_needed = {
        "IR": ["IR", "IS", "IM"],
        "IIR": ["IIR", "IIM"],
        "IIIR": ["IIIR", "IIIM"]
    }
    if system_type in types_needed:
        required_types = types_needed[system_type]
        present_types = [row["type"] for row in components]

        # Ensure all required types are present
        result = all(r_type in present_types for r_type in required_types)
        return result
    else:
        # For IV and IIG, no additional types are needed
        return True


def format_results(results):
    formatted = []
    for system in results:
        formatted.append("-------------------------------------------------------------------------")
        for component in system:
            formatted.append(
                f"{component['type']}\t{component['module']}\t{component['genome.acc']}\t{component['protein.id']}"
            )
        formatted.append("-------------------------------------------------------------------------")
    return formatted


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py neighborhood.txt")
        sys.exit(1)

    filename = sys.argv[1]
    data = parse_file(filename)
    defense_systems = find_defense_systems(data)
    results = format_results(defense_systems)

    with open("results.tsv", "w") as output_file:
        for line in results:
            output_file.write(line + "\n")
