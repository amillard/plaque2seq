import os
import shutil
import fnmatch
import argparse

# Directories with literal names to delete
TARGET_DIRS = {
    'taxmyphage', 'pipeline_state' ,'clean', 'bam', 'polish1', 'polish2', 'tmp', 'misc', 'unicycler*' ,  'K21', 'K33', 'K55', 'K77', 'K99', 'K127',  
}

# Directory name patterns to match (e.g., 'flye*')
TARGET_DIR_PATTERNS = ['flye*']

# File patterns to delete (glob style)
DELETE_PATTERNS = [
    'SQK-RBK114-96_barcode*.fa',
    'SQK-RBK114-96_barcode*.fa.split.fasta',
    'SQK-RBK114-96_barcode*.fa.non_chimera.fasta',
    'SQK-RBK114-24_barcode*.fa',
    'SQK-RBK114-24_barcode*.fa.split.fasta',
    'SQK-RBK114-24_barcode*.fa.non_chimera.fasta',
    'scaffolds.fasta',
    'scaffolds.paths',
    'contigs.paths',
    '*yaml',
    'assembly_graph_with_scaffolds.gfa',
    'before_rr.fasta',
]

# Suffixes of files to keep
KEEP_SUFFIXES = ('.output.fasta', '.output.fasta.gz')

# Files to delete only inside 'checkv' directories
CHECKV_TARGET_FILES = {'proviruses.fna', 'viruses.fna'}

def clean_directory(main_directory, dry_run=False):
    for dirpath, dirnames, filenames in os.walk(main_directory, topdown=True):
        current_dirname = os.path.basename(dirpath)

        # Remove matching subdirectories
        for dirname in dirnames[:]:
            full_path = os.path.join(dirpath, dirname)

            if dirname in TARGET_DIRS or any(fnmatch.fnmatch(dirname, pat) for pat in TARGET_DIR_PATTERNS):
                _delete(full_path, dry_run, is_dir=True)
                dirnames.remove(dirname)

        # File-level deletion
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)

            # Check if it's in a 'checkv' directory and matches target name
            if current_dirname == 'checkv' and filename in CHECKV_TARGET_FILES:
                _delete(file_path, dry_run)
                continue

            if filename.endswith(KEEP_SUFFIXES):
                continue

            if any(fnmatch.fnmatch(filename, pattern) for pattern in DELETE_PATTERNS):
                _delete(file_path, dry_run)

    print("\n✅ DRY RUN complete. No changes were made." if dry_run else "\n✅ Cleanup complete.")

def _delete(path, dry_run, is_dir=False):
    if dry_run:
        print(f"[DRY RUN] Would delete {'directory' if is_dir else 'file'}: {path}")
    else:
        print(f"Deleting {'directory' if is_dir else 'file'}: {path}")
        if is_dir:
            shutil.rmtree(path)
        else:
            os.remove(path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean up target directories and sequencing files.")
    parser.add_argument("directory", help="Main directory to start cleanup from")
    parser.add_argument("--dry-run", action="store_true", help="Preview what would be deleted without actually deleting anything")

    args = parser.parse_args()
    clean_directory(args.directory, dry_run=args.dry_run)
