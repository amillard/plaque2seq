#!/usr/bin/env python3

# splits long reads into short reads to allow future processing
# splits into 300 bp fragments by default
# KEEPS fragments shorter than chunk size (e.g. tail fragments)

import os
import glob
import gzip
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_barcode_mapping(mapping_file: str) -> dict:
    mapping = {}
    with open(mapping_file) as fh:
        next(fh)
        for line in fh:
            barcode, sample = line.strip().split()
            mapping[barcode] = sample
    return mapping


def _test_gzip(path: str) -> bool:
    try:
        subprocess.run(
            ["gzip", "-t", path],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def chunk_file(input_file: str, sample_name: str, chunk_size: int = 300, output_dir: str = ".") -> str:
    if input_file.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
        fmt, ext = "fastq", "fq.gz"
    elif input_file.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
        fmt, ext = "fasta", "fa.gz"
    else:
        raise ValueError(f"Unrecognized extension: {input_file!r}")

    os.makedirs(output_dir, exist_ok=True)
    out_name = os.path.join(output_dir, f"{sample_name}_short_{chunk_size}.{ext}")

    if os.path.exists(out_name):
        if _test_gzip(out_name):
            print(f"Skipping {input_file!r}: {out_name!r} exists & is valid.")
            return out_name
        else:
            print(f"Corrupted archive {out_name!r}, removing & retrying.")
            os.remove(out_name)

    for attempt in (1, 2):
        with gzip.open(out_name, "wt") as out_fh:
            with (gzip.open(input_file, "rt") if input_file.endswith(".gz") else open(input_file, "rt")) as in_fh:
                for rec in SeqIO.parse(in_fh, fmt):
                    seq = rec.seq
                    quals = rec.letter_annotations.get("phred_quality", None)

                    for i in range(0, len(seq), chunk_size):
                        chunk = seq[i : i + chunk_size]

                        # --- CHANGED LOGIC STARTS HERE ---
                        if fmt == "fastq":
                            if not quals:
                                continue
                            chunk_quals = quals[i : i + len(chunk)]
                            if len(chunk_quals) != len(chunk):
                                continue
                        # --- CHANGED LOGIC ENDS HERE ---

                        new_id = f"{rec.id}_chunk{i}"
                        new_rec = SeqRecord(chunk, id=new_id, description="")

                        if fmt == "fastq":
                            new_rec.letter_annotations["phred_quality"] = chunk_quals

                        SeqIO.write(new_rec, out_fh, fmt)

        if _test_gzip(out_name):
            print(f" Split into fragments {input_file!r} → {out_name!r}")
            return out_name
        else:
            print(f"Attempt {attempt} failed gzip test for {out_name!r}")
            os.remove(out_name)

    raise IOError(f"Failed to write valid gzip for {out_name!r} after retrying")


# This function runs in a separate process for each file
def worker_wrapper(infile, mapping, chunk_size, output_dir):
    basename = os.path.basename(infile)
    parts = basename.split("_")
    if len(parts) < 2 or not parts[1].startswith("barcode"):
        raise ValueError(f"Unexpected filename format: {basename}")
    bc = parts[1].split(".")[0]  # e.g., barcode34
    sample = mapping.get(bc)
    if not sample:
        raise ValueError(f"No sample mapping for barcode {bc}")
    return chunk_file(infile, sample, chunk_size, output_dir)


# assumes the standard format of SQK
def process_directory(
    directory: str,
    mapping_file: str,
    chunk_size: int = 300,
    pattern: str = "SQK-RBK114-96_*.fastq.gz",
    workers: int = None,
    output_dir: str = ".",
) -> list:
    mapping = read_barcode_mapping(mapping_file)
    files = glob.glob(os.path.join(directory, pattern))
    results = []

    with ProcessPoolExecutor(max_workers=workers) as exe:
        futures = {
            exe.submit(worker_wrapper, f, mapping, chunk_size, output_dir): f for f in files
        }

        for future in as_completed(futures):
            infile = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Failed to process {infile}: {e}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Chunk long FASTQ/FASTA files into fixed-length fragments.")
    parser.add_argument("-d", "--directory", required=True, help="Input directory with FASTQ/FASTA files.")
    parser.add_argument("-m", "--mapping", required=True, help="Barcode-to-sample mapping file.")
    parser.add_argument("-c", "--chunk_size", type=int, default=300, help="Chunk size (default: 300).")
    parser.add_argument("-p", "--pattern", default="SQK-RBK114-96_*.fastq.gz", help="Glob pattern to match files.")
    parser.add_argument("-w", "--workers", type=int, default=None, help="Number of parallel workers.")
    parser.add_argument("-o", "--output", default=".", help="Output directory.")

    args = parser.parse_args()

    process_directory(
        directory=args.directory,
        mapping_file=args.mapping,
        chunk_size=args.chunk_size,
        pattern=args.pattern,
        workers=args.workers,
        output_dir=args.output,
    )


if __name__ == "__main__":
    main()
