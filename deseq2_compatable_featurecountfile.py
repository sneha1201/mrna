#!/usr/bin/env python3

import argparse
import pandas as pd
import re

def clean_sample_name(name):
    name = re.sub(r'^.*/', '', name)  # remove path
    name = re.sub(r'\.bam$|\.sam$|\.txt$|\.tsv$|\.csv$', '', name)  # remove extensions
    name = re.sub(r'_sorted$', '', name)  # remove trailing _sorted
    return name

def main():
    parser = argparse.ArgumentParser(description="Clean featureCounts file for DESeq2")
    parser.add_argument("-i", "--input", required=True, help="Input featureCounts file")
    parser.add_argument("-o", "--output", required=True, help="Output cleaned file")
    args = parser.parse_args()

    with open(args.input, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith("#") and line.strip() != ""]

    # Find header by scanning first non-comment line with 'Geneid'
    header_line = next((line for line in lines if 'Geneid' in line), None)
    if header_line is None:
        raise ValueError("❌ Could not find header line containing 'Geneid'.")

    header = header_line.strip().split("\t")
    print(f"🔍 Header column count: {len(header)}")

    header_index = lines.index(header_line)
    data_lines = lines[header_index + 1:]

    matched_rows = []
    mismatched_count = 0

    for i, line in enumerate(data_lines, start=header_index + 2):
        cols = line.split("\t")
        if len(cols) == len(header):
            matched_rows.append(cols)
        else:
            mismatched_count += 1
            print(f"⚠️ Line {i}: {len(cols)} columns (expected {len(header)})")

    if not matched_rows:
        raise ValueError("❌ No valid data rows matched header length.")

    print(f"✅ Rows matched: {len(matched_rows)}")
    print(f"⚠️ Rows skipped due to mismatch: {mismatched_count}")

    df = pd.DataFrame(matched_rows, columns=header)

    if "Geneid" not in df.columns:
        raise ValueError("❌ 'Geneid' column missing in cleaned header.")

    # Keep only relevant columns
    df = df.loc[:, ["Geneid"] + df.columns[6:].tolist()]
    df.columns = ["Geneid"] + [clean_sample_name(col) for col in df.columns[1:]]

    df.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Cleaned file saved to: {args.output}")

if __name__ == "__main__":
    main()

