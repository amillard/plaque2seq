print("Will now create mash database(s)")
from Bio import SeqIO

def _fasta_has_records(fp: str) -> bool:
    try:
        if not (os.path.isfile(fp) and os.path.getsize(fp) > 0):
            return False
        with open(fp) as fh:
            return next(SeqIO.parse(fh, "fasta"), None) is not None
    except Exception:
        return False

try:
    candidates = []
    for fp in (path_genomes, path_genomes_exclude_refseq):
        if _fasta_has_records(fp):
            candidates.append(fp)
        else:
            print(f"[mash] Skipping (missing or no records): {fp}")

    if not candidates:
        print("[mash] Nothing to sketch; skipping Mash.")
    else:
        for fp in candidates:
            mymodule.create_mash(fp, mash_path, threads=args.num_cpus, dry_run=False)
except Exception as e:
    print(f"Mash creation failed: {e}")
