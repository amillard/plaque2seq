# Output file for summary

import os 

import argparse
parser = argparse.ArgumentParser(description="Combine all check files into a single file .")
parser.add_argument("--output_file", default="combined_quality_summary.tsv",
                        help="Path to log file (default: checkv_processing.log)")
args = parser.parse_args()


#output_file = "combined_quality_summary.tsv"

output_file = args.output_file

print (f"will concatenate all checkv files into a single file {output_file} for further processing ")


with open(output_file, "w") as out_f:
    out_f.write("Directory\tOriginal_Line\n")  # Header line

    for root, dirs, files in os.walk("."):  # Start from current directory
        if 'checkv' in dirs:
            checkv_path = os.path.join(root, 'checkv')
            quality_file = os.path.join(checkv_path, 'quality_summary.tsv')
            relative_dir = os.path.relpath(root)

            if os.path.isfile(quality_file):
                found = False
                with open(quality_file, "r") as qf:
                    for line in qf:
                        if 'High' in line or 'Medium' in line:
                            out_f.write(f"{relative_dir}\t{line.strip()}\n")
                            found = True
                if not found:
                    out_f.write(f"{relative_dir}\tCheck this directory out\n")
            else:
                out_f.write(f"{relative_dir}\tNo checkv produced\n")


print (f"complete concatenation of all checkv quality summary files ")
