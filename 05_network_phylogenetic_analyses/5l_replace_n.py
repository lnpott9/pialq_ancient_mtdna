#!/usr/bin/env python3

# =============================================================================
# Title: 5l_replace_n.py
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Python script to create map showing geographic metadata for reference database of mtDNA sequences used for networks
# Usage: python 5l_replace_n.py <input> <output>
# =============================================================================

import argparse

def replace_N_with_missing(infile, outfile):
    with open(infile, "r") as f:
        lines = f.readlines()

    new_lines = []
    in_matrix = False

    for line in lines:
        stripped = line.strip()

        # detect start/end of the MATRIX section
        if stripped.upper().startswith("MATRIX"):
            in_matrix = True
            new_lines.append(line)
            continue
        if in_matrix and stripped.upper().startswith(";"):
            in_matrix = False
            new_lines.append(line)
            continue

        if in_matrix and stripped:
            parts = stripped.split()
            if len(parts) >= 2:
                name, seq = parts[0], parts[1]
                # Replace Ns (upper/lower) with ?
                seq = seq.replace("N", "?").replace("n", "?")
                new_line = f"{name}\t{seq}\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Write out modified nexus file
    with open(outfile, "w") as f:
        f.writelines(new_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Replace Ns in NEXUS DNA sequences with ?"
    )
    parser.add_argument("infile", help="Input NEXUS file")
    parser.add_argument("outfile", help="Output NEXUS file with Ns replaced by ?")

    args = parser.parse_args()
    replace_N_with_missing(args.infile, args.outfile)


if __name__ == "__main__":
    main()