#!/usr/bin/env python3
"""
Build an nf-core/bacass samplesheet from folders of FASTQ files.

bacass samplesheet columns (6): ID, R1, R2, LongFastQ, Fast5, GenomeSize
  - ID         : unique; must start with a letter; letters/numbers/underscores only
  - R1, R2     : short (Illumina) paired-end reads
  - LongFastQ  : long reads (Nanopore/PacBio)
  - Fast5      : directory of FAST5 files (left as NA here)
  - GenomeSize : a number ending in 'm', e.g. 4.6m
Empty cells are written as NA. A directory is scanned NON-recursively (only files
directly inside it - subdirectories are not searched); a quoted glob such as
'short/*.fastq.gz' also works. Local paths or gs:// prefixes (see GOOGLE CLOUD
STORAGE below for gs:// auth checks).

==============================================================================
ARGUMENTS
==============================================================================
Required (every run):
  fastq_location        Positional. A directory (scanned non-recursively) or a
                        quoted glob like 'short/*.fastq.gz', local or gs://. Holds
                        the SHORT reads by default, or the LONG reads with --long-only.
  output                Positional. Path to write the sheet. The extension sets
                        the delimiter: .csv -> comma, .tsv -> tab.

Optional (and when to use them):
  --long-only           Use ONLY when you have long reads and no short reads
                        (ONT/PacBio-only assembly). Makes fastq_location the
                        long reads. Cannot be combined with --long-fastq-dir.
  --long-fastq-dir DIR  Use ONLY for hybrid runs (short + long). DIR is the dir
                        or gs:// prefix holding the long reads; they are matched
                        to the short samples by name. Cannot be combined with
                        --long-only.
  --genome-size SIZE    Optional always. Recommended for --long-only and hybrid
                        runs (long-read assemblers/QC want the estimate); unused
                        by Illumina-only short-read assembly. Format: number
                        ending in 'm' (e.g. 4.6m). Default: NA.
  --include GLOB        Optional. Keep only files whose basename matches GLOB
                        (e.g. 'SampleA*') to pull one pair (or a subset) out of a
                        location holding many. Repeatable. Quote it so the shell
                        does not expand it first.

Mode  ->  required args                              ->  bacass launch flag
  Illumina-only   fastq_location, output                  --assembly_type short
  Long-only       fastq_location, output, --long-only     --assembly_type long
  Hybrid          fastq_location, output, --long-fastq-dir --assembly_type hybrid
==============================================================================

------------------------------------------------------------------------------
GOOGLE CLOUD STORAGE (gs://) INPUTS
------------------------------------------------------------------------------
A gs:// input is listed with `gsutil`, so the Google Cloud SDK must be installed
and authenticated on the machine RUNNING this script (your own identity, separate
from the pipeline's service account) with read access - roles/storage.objectViewer
- to the bucket.

Confirm you are authenticated to the right project and can see the bucket BEFORE
running the script:

    gcloud auth list                      # which account is active (* = active)
    gcloud config get-value project       # which project is currently set
    gcloud auth login                     # sign in if no account is active
    gcloud config set project PROJECT_ID  # switch to the project that owns the bucket
    gsutil ls gs://YOUR_BUCKET/           # list the bucket - the real read test

If that last command lists objects, the script can read that location. AccessDenied
or 401 means fix auth/permissions first; gsutil "not found" means install the SDK.

Usage:
    # Illumina-only -> CSV (genome size not needed for short-read assembly)
    ./make_bacass_samplesheet.py gs://my-bucket/short/ samplesheet.csv

    # Hybrid: short reads here + long reads elsewhere, with genome size
    ./make_bacass_samplesheet.py gs://my-bucket/short/ samplesheet.csv \
        --long-fastq-dir gs://my-bucket/long/ --genome-size 4.6m

    # ONT/PacBio-only (point the main location at the long reads)
    ./make_bacass_samplesheet.py gs://my-bucket/long/ samplesheet.tsv \
        --long-only --genome-size 4.6m

    # Point at a glob of specific files instead of a whole directory
    ./make_bacass_samplesheet.py 'short/*.fastq.gz' samplesheet.csv

    # One pair out of a location holding many (glob on file basename)
    ./make_bacass_samplesheet.py gs://my-bucket/short/ sampleA.csv --include 'SampleA*'
"""
import argparse
import csv
import fnmatch
import glob
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# Short-read pairing: SAMPLE_R1_001.fastq.gz, SAMPLE_R1.fastq.gz, SAMPLE_1.fastq.gz
PAIR_RE = re.compile(r"^(?P<sample>.+?)_R?(?P<read>[12])(?:_\d{3})?\.(?:fastq|fq)\.gz$")
# Tokens stripped from long-read filenames to recover the sample key
LONG_TOKENS = ("_longreads", "_nanopore", "_pacbio", "_minion", "_long", "_ont", "_lr")
FASTQ_EXTS = (".fastq.gz", ".fq.gz")
ID_VALID = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def run_gsutil_ls(target):
    """Run `gsutil ls target`; return stdout lines, exiting cleanly on failure."""
    try:
        result = subprocess.run(["gsutil", "ls", target],
                                capture_output=True, text=True)
    except FileNotFoundError:
        sys.exit("gsutil not found on PATH; install and authenticate the "
                 "Google Cloud SDK (gcloud auth login) to read gs:// paths")
    if result.returncode != 0:
        if "matched no objects" in result.stderr:
            return []  # caller reports the friendly "No FASTQ files found"
        sys.exit(f"gsutil ls failed for {target}:\n{result.stderr.strip()}")
    return result.stdout.splitlines()


def resolve_inputs(location, patterns=None):
    """FASTQ paths/URIs from a directory (NON-recursive) or a glob pattern.

    `location` is either a directory - only files directly inside it are used,
    subdirectories are not searched - or a shell glob like 'short/*.fastq.gz'.
    `patterns` is the optional --include filter, matched against basenames.
    """
    is_glob = any(ch in location for ch in "*?[")
    if location.startswith("gs://"):
        # Directory: trailing slash lists immediate children only (no recursion).
        # Glob: hand the pattern to gsutil as-is.
        target = location if is_glob else location.rstrip("/") + "/"
        files = [ln.strip() for ln in run_gsutil_ls(target)
                 if ln.strip().endswith(FASTQ_EXTS)]
    elif is_glob:
        files = [f for f in glob.glob(location) if f.endswith(FASTQ_EXTS)]
    else:
        files = [str(p) for p in Path(location).glob("*")  # immediate entries only
                 if p.name.endswith(FASTQ_EXTS)]
    if patterns:
        files = [f for f in files
                 if any(fnmatch.fnmatch(f.split("/")[-1], pat) for pat in patterns)]
    return files


def pair_short_reads(files):
    """Group paired/single short reads by sample name into r1/r2 slots."""
    samples = defaultdict(dict)
    for uri in files:
        name = uri.split("/")[-1]
        match = PAIR_RE.match(name)
        if not match:
            print(f"WARNING: '{name}' has no R1/R2 read number; skipping as a short read",
                  file=sys.stderr)
            continue
        samples[match.group("sample")][f"r{match.group('read')}"] = uri
    return samples


def long_sample_key(name):
    """Recover a sample key from a long-read filename (strip ext + long-read tags)."""
    base = name
    for ext in FASTQ_EXTS:
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    for token in LONG_TOKENS:
        if base.lower().endswith(token):
            base = base[: -len(token)]
            break
    return base


def collect_long_reads(files):
    """Map sample key -> long-read URI."""
    longs = {}
    for uri in files:
        key = long_sample_key(uri.split("/")[-1])
        if key in longs:
            print(f"WARNING: multiple long-read files map to '{key}'; keeping the first",
                  file=sys.stderr)
            continue
        longs[key] = uri
    return longs


def sanitize_id(raw):
    """Return a bacass-valid ID (letter start; letters/numbers/underscores only)."""
    if ID_VALID.match(raw):
        return raw
    fixed = re.sub(r"[^A-Za-z0-9_]", "_", raw)
    if not re.match(r"^[A-Za-z]", fixed):
        fixed = f"sample_{fixed}"
    print(f"WARNING: sample ID '{raw}' is not bacass-valid; renamed to '{fixed}'",
          file=sys.stderr)
    return fixed


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # --- Required positionals (supplied on every run) ---
    parser.add_argument(
        "fastq_location",
        help="REQUIRED. Dir or gs:// prefix to scan. Short reads by default; "
             "long reads when --long-only is set.")
    parser.add_argument(
        "output",
        help="REQUIRED. Output samplesheet path; .csv -> comma, .tsv -> tab.")

    # --- Read-type selection (OPTIONAL; choose at most one) ---
    #   neither flag         -> Illumina short-read-only sheet  (--assembly_type short)
    #   --long-only          -> long-read-only sheet            (--assembly_type long)
    #   --long-fastq-dir DIR -> hybrid sheet                    (--assembly_type hybrid)
    # The two flags are alternatives, so they share a mutually exclusive group:
    # passing both is rejected by argparse rather than silently resolved.
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--long-only", action="store_true",
        help="OPTIONAL. Use when there are ONLY long reads (ONT/PacBio). Treats "
             "fastq_location as the long reads. Not combinable with --long-fastq-dir.")
    mode.add_argument(
        "--long-fastq-dir", default=None, metavar="DIR",
        help="OPTIONAL. Use for HYBRID runs. Dir or gs:// prefix holding the long "
             "reads, matched to short samples by name. Not combinable with --long-only.")

    # --- Other optional settings ---
    parser.add_argument(
        "--genome-size", default="NA", metavar="SIZE",
        help="OPTIONAL. Recommended for --long-only and hybrid runs; unused by "
             "Illumina-only assembly. Number ending in 'm' (e.g. 4.6m). Default: NA.")
    parser.add_argument(
        "--include", action="append", metavar="GLOB",
        help="OPTIONAL. Keep only files whose basename matches this glob (e.g. "
             "'SampleA*') to pick one pair/subset from a crowded location. "
             "Repeatable. Quote it.")

    args = parser.parse_args()

    # GenomeSize is free-form text in the sheet; nudge toward bacass's 'm' notation.
    if args.genome_size != "NA" and not args.genome_size.strip().endswith("m"):
        print(f"WARNING: --genome-size '{args.genome_size}' should end in 'm' "
              "(e.g. 4.6m) for bacass", file=sys.stderr)

    rows = {}  # sanitized_id -> {R1, R2, LongFastQ}

    if args.long_only:
        # Long-read-only mode: fastq_location IS the long reads; R1/R2 stay NA.
        long_files = resolve_inputs(args.fastq_location, args.include)
        if not long_files:
            sys.exit(f"No FASTQ files found under {args.fastq_location}")
        for key, uri in collect_long_reads(long_files).items():
            rows[sanitize_id(key)] = {"R1": "NA", "R2": "NA", "LongFastQ": uri}
    else:
        # Short-read mode, optionally upgraded to hybrid by --long-fastq-dir.
        short_files = resolve_inputs(args.fastq_location, args.include)
        if not short_files:
            sys.exit(f"No FASTQ files found under {args.fastq_location}")
        for sample, reads in pair_short_reads(short_files).items():
            rows[sanitize_id(sample)] = {
                "R1": reads.get("r1", "NA"),
                "R2": reads.get("r2", "NA"),
                "LongFastQ": "NA",
            }
        if args.long_fastq_dir:
            # Hybrid: attach long reads to matching short samples (or add new rows).
            for key, uri in collect_long_reads(
                    resolve_inputs(args.long_fastq_dir, args.include)).items():
                sid = sanitize_id(key)
                if sid in rows:
                    rows[sid]["LongFastQ"] = uri
                else:
                    rows[sid] = {"R1": "NA", "R2": "NA", "LongFastQ": uri}

    if not rows:
        sys.exit("No samples assembled; check your inputs")

    delimiter = "\t" if args.output.endswith(".tsv") else ","
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        writer.writerow(["ID", "R1", "R2", "LongFastQ", "Fast5", "GenomeSize"])
        for sample_id in sorted(rows):
            r = rows[sample_id]
            writer.writerow([sample_id, r["R1"], r["R2"], r["LongFastQ"],
                             "NA", args.genome_size])

    print(f"Wrote {len(rows)} samples to {args.output}")


if __name__ == "__main__":
    main()
