#!/usr/bin/env python3

import os
import csv

def parse_checkv_quality_summary(file_path):
    """
    Parse checkv/quality_summary.tsv and count number of contigs with
    High-quality, Medium-quality, or Complete status in the 'checkv_quality' column.
    """
    counts = {"High-quality": 0, "Medium-quality": 0, "Complete": 0}
    with open(file_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            quality = row.get("checkv_quality", "")
            if "High" in quality:
                counts["High-quality"] += 1
            if "Medium" in quality:
                counts["Medium-quality"] += 1
            if "Complete" in quality:
                counts["Complete"] += 1
    return counts

def find_quality_summaries(base_dir="."):
    """
    Walk the directory tree to find all checkv/quality_summary.tsv files.
    """
    found = []
    for root, dirs, files in os.walk(base_dir):
        if "checkv" in dirs:
            quality_file = os.path.join(root, "checkv", "quality_summary.tsv")
            if os.path.isfile(quality_file):
                found.append((os.path.relpath(root), quality_file))
    return found

def main():
    output_file = "checkv_quality_summary_counts.tsv"
    summaries = find_quality_summaries(".")

    with open(output_file, "w", newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(["Sample_Directory", "Complete",   "High-quality", "Medium-quality"])

        for sample_dir, quality_path in summaries:
            try:
                counts = parse_checkv_quality_summary(quality_path)
                writer.writerow([
                    sample_dir,
                    counts["Complete"],
                    counts["High-quality"],
                    counts["Medium-quality"],
                ])
            except Exception as e:
                print(f"Failed to process {quality_path}: {e}")

    print(f"Summary written to {output_file}")

if __name__ == "__main__":
    main()
