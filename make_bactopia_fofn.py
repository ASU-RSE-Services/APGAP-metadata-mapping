#!/usr/bin/env python3
"""
Build a Bactopia FOFN (file of filenames / samplesheet) from folders of reads.

The FOFN is a TAB-delimited table with 7 columns, unchanged from Bactopia v3.0.0
onward (the v4 release changed how Bactopia runs - Nextflow 26.04 / nf-bactopia -
not the input schema), so this one script serves both v3.x and v4.x:

    sample  runtype  genome_size  species  r1  r2  extra

  - sample      : unique name used to label all outputs (Bactopia requires unique)
  - runtype     : paired-end | single-end | ont | hybrid | short_polish
  - genome_size : integer bp (0 = Bactopia's default handling; or pass a value)
  - species     : species name (UNKNOWN_SPECIES if not given)
  - r1          : R1 (paired-end), the single-end read, or the ONT read (runtype ont)
  - r2          : R2 for paired-end / hybrid / short_polish; empty otherwise
  - extra       : the long read for hybrid / short_polish; empty otherwise

Empty cells are left blank (tab-separated), matching `bactopia prepare`. Works on
local directories or gs:// prefixes (see GOOGLE CLOUD STORAGE below for gs:// auth
checks). Output is always tab-delimited with clean LF line endings.

NOTE: verify the header on your install with `bactopia prepare --help` if you want
to be certain v4 matches; if it differs, this becomes a one-line `--bactopia-version`
branch rather than a second script.

==============================================================================
ARGUMENTS
==============================================================================
Required (every run):
  fastq_location        Positional. A directory (scanned non-recursively - only
                        files directly inside it) or a quoted glob like
                        'short/*.fastq.gz', local or gs://. Holds the SHORT reads
                        by default, or the ONT reads with --ont.
  output                Positional. Path to write the FOFN (tab-delimited; .tsv
                        or .txt by convention).

Optional (and when to use them):
  --ont                 Use when the reads are Nanopore/PacBio long reads with NO
                        short reads (runtype 'ont'). Makes fastq_location the long
                        reads. Mutually exclusive with --long-fastq-dir.
  --long-fastq-dir DIR  Use for hybrid runs (short + long). DIR is the dir or gs://
                        prefix holding the long reads, matched to short samples by
                        name. Produces runtype 'hybrid'. Mutually exclusive with
                        --ont.
  --short-polish        Only with --long-fastq-dir. Emit runtype 'short_polish'
                        (ONT primary, Illumina used to polish) instead of 'hybrid'.
  --genome-size SIZE    Integer bp applied to every sample. Default: 0
                        (Bactopia's default genome-size handling).
  --species NAME        Species applied to every sample. Default: UNKNOWN_SPECIES.
  --include GLOB        Keep only files whose basename matches GLOB (e.g.
                        'SampleA*') to pull one sample/subset from a crowded
                        location. Repeatable. Quote it.

Mode -> runtype -> Bactopia launch
  short paired-end   paired-end    bactopia --samples <fofn> ...
  short single-end   single-end    bactopia --samples <fofn> ...
  --ont              ont           bactopia --samples <fofn> ...
  --long-fastq-dir   hybrid        bactopia --samples <fofn> ...
  + --short-polish   short_polish  bactopia --samples <fofn> ...
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
    # short-read paired-end -> TSV
    ./make_bactopia_fofn.py gs://my-bucket/short/ samples.tsv \
        --genome-size 2800000 --species "Staphylococcus aureus"

    # ONT-only
    ./make_bactopia_fofn.py gs://my-bucket/ont/ samples.tsv --ont

    # hybrid (short here + long elsewhere)
    ./make_bactopia_fofn.py gs://my-bucket/short/ samples.tsv \
        --long-fastq-dir gs://my-bucket/long/

    # one sample out of many
    ./make_bactopia_fofn.py gs://my-bucket/short/ one.tsv --include 'SampleA*'

    # a glob of specific files instead of a whole directory
    ./make_bactopia_fofn.py 'short/*.fastq.gz' samples.tsv
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

FOFN_HEADER = ["sample", "runtype", "genome_size", "species", "r1", "r2", "extra"]
PAIR_RE = re.compile(r"^(?P<sample>.+?)_R?(?P<read>[12])(?:_\d{3})?\.(?:fastq|fq)\.gz$")
LONG_TOKENS = ("_longreads", "_nanopore", "_pacbio", "_minion", "_long", "_ont", "_lr")
FASTQ_EXTS = (".fastq.gz", ".fq.gz")


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
            return []
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
    """Group paired/single short reads by sample into r1/r2 slots."""
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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # --- Required positionals (every run) ---
    parser.add_argument(
        "fastq_location",
        help="REQUIRED. Dir or gs:// prefix to scan. Short reads by default; "
             "long reads when --ont is set.")
    parser.add_argument(
        "output", help="REQUIRED. Output FOFN path (tab-delimited).")

    # --- Read-type selection (OPTIONAL; choose at most one) ---
    #   neither flag          -> short reads -> runtype paired-end / single-end
    #   --ont                 -> long-read-only -> runtype ont
    #   --long-fastq-dir DIR  -> hybrid -> runtype hybrid (or short_polish)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--ont", action="store_true",
        help="OPTIONAL. Reads are Nanopore/PacBio long reads, no short reads "
             "(runtype 'ont'). Not combinable with --long-fastq-dir.")
    mode.add_argument(
        "--long-fastq-dir", default=None, metavar="DIR",
        help="OPTIONAL. Dir or gs:// prefix of long reads for a HYBRID run "
             "(runtype 'hybrid'), matched to short samples by name. Not combinable "
             "with --ont.")

    # --- Other optional settings ---
    parser.add_argument(
        "--short-polish", action="store_true",
        help="OPTIONAL. With --long-fastq-dir, emit runtype 'short_polish' instead "
             "of 'hybrid' (ONT primary, Illumina polishes).")
    parser.add_argument(
        "--genome-size", default="0", metavar="SIZE",
        help="OPTIONAL. Integer bp for every sample. Default: 0.")
    parser.add_argument(
        "--species", default="UNKNOWN_SPECIES", metavar="NAME",
        help="OPTIONAL. Species for every sample. Default: UNKNOWN_SPECIES.")
    parser.add_argument(
        "--include", action="append", metavar="GLOB",
        help="OPTIONAL. Keep only files whose basename matches this glob (e.g. "
             "'SampleA*'). Repeatable. Quote it.")

    args = parser.parse_args()

    if args.short_polish and not args.long_fastq_dir:
        sys.exit("--short-polish requires --long-fastq-dir")

    gsize, species = args.genome_size, args.species
    rows = {}  # sample -> [runtype, r1, r2, extra]

    def add_row(sample, runtype, r1, r2, extra):
        if sample in rows:
            print(f"WARNING: duplicate sample '{sample}'; keeping the first",
                  file=sys.stderr)
            return
        rows[sample] = [runtype, r1, r2, extra]

    if args.ont:
        # ONT-only: each long read is its own sample, runtype 'ont', read in r1.
        ont_files = resolve_inputs(args.fastq_location, args.include)
        if not ont_files:
            sys.exit(f"No FASTQ files found under {args.fastq_location}")
        for uri in ont_files:
            add_row(long_sample_key(uri.split("/")[-1]), "ont", uri, "", "")
    else:
        # Short reads, optionally upgraded to hybrid/short_polish via --long-fastq-dir.
        short_files = resolve_inputs(args.fastq_location, args.include)
        if not short_files:
            sys.exit(f"No FASTQ files found under {args.fastq_location}")
        longs = {}
        if args.long_fastq_dir:
            longs = collect_long_reads(resolve_inputs(args.long_fastq_dir, args.include))
        hybrid_runtype = "short_polish" if args.short_polish else "hybrid"

        for sample, reads in pair_short_reads(short_files).items():
            r1 = reads.get("r1", "")
            r2 = reads.get("r2", "")
            if not r1:
                print(f"WARNING: sample '{sample}' has R2 but no R1; skipping",
                      file=sys.stderr)
                continue
            long_uri = longs.pop(sample, None)
            if args.long_fastq_dir and long_uri:
                add_row(sample, hybrid_runtype, r1, r2, long_uri)
            elif r2:
                add_row(sample, "paired-end", r1, r2, "")
            else:
                add_row(sample, "single-end", r1, "", "")

        # Long reads with no matching short sample (hybrid mode): warn, they are dropped.
        for key in longs:
            print(f"WARNING: long read for '{key}' has no matching short-read sample; "
                  "skipped", file=sys.stderr)

    if not rows:
        sys.exit("No samples assembled; check your inputs")

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FOFN_HEADER)
        for sample in sorted(rows):
            runtype, r1, r2, extra = rows[sample]
            writer.writerow([sample, runtype, gsize, species, r1, r2, extra])

    print(f"Wrote {len(rows)} samples to {args.output}")


if __name__ == "__main__":
    main()
