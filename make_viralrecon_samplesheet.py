#!/usr/bin/env python3
"""
Build an nf-core/viralrecon samplesheet for Illumina OR Nanopore data.

viralrecon runs ONE platform per execution, and each platform needs a different
samplesheet. This script writes whichever one you ask for:

  --platform illumina   (DEFAULT)  ->  columns: sample,fastq_1,fastq_2
  --platform nanopore              ->  columns: sample,barcode

Illumina and Nanopore CANNOT be combined in one viralrecon run, so this script
produces one platform's sheet per call. Output is always CSV (viralrecon requires
comma-separated input).

------------------------------------------------------------------------------
INPUT (first positional argument)
------------------------------------------------------------------------------
Every directory or glob below may be a LOCAL path OR a gs:// Google Cloud Storage
URI. gs:// inputs are listed with `gsutil`, and the gs:// paths are written into
the samplesheet as-is (what viralrecon needs on Google Cloud Batch). See GOOGLE
CLOUD STORAGE below for the auth checks to run first.

illumina:
  - a DIRECTORY of FASTQ files, e.g.  reads/  or  gs://my-bucket/reads/
    Only files directly inside it are used; SUBDIRECTORIES ARE NOT SEARCHED.
  - or a GLOB, e.g.  'reads/*.fastq.gz'  or  'gs://my-bucket/reads/*.fastq.gz'
    Quote it, otherwise your shell expands it before the script sees it.
  Reads must be gzipped (.fastq.gz / .fq.gz). R1/R2 are paired by filename;
  a sample with only R1 becomes single-end (empty fastq_2).

nanopore:
  - the FASTQ_DIR, e.g.  fastq_pass/  or  gs://my-bucket/fastq_pass/
    A directory whose immediate subdirectories are barcodeNN/. The script lists
    those barcode folders; it does not read the FASTQs inside them (viralrecon
    locates the reads itself at run time via --fastq_dir). Pass the directory,
    not a glob.

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

------------------------------------------------------------------------------
ARGUMENTS
------------------------------------------------------------------------------
Required (every run):
  input              Directory or, for illumina, a glob (see INPUT above).
  output             Output samplesheet path (written as CSV).

Optional:
  --platform NAME    'illumina' (DEFAULT) or 'nanopore'. Omit it for Illumina;
                     you only need it for Nanopore data.
  --include GLOB     illumina only. Keep only files whose basename matches GLOB
                     (e.g. 'SampleA*'); repeatable. Another way to narrow things
                     when you pass a whole directory. Rejected with nanopore.
  --barcode-map FILE nanopore only. 2-column CSV/TSV 'barcode,sample' giving a
                     real name per barcode; barcodes not listed keep their folder
                     name (e.g. barcode01). Rejected with illumina.

------------------------------------------------------------------------------
WHAT YOU RUN AFTER THIS SCRIPT (reference only - this script does not run it)
------------------------------------------------------------------------------
This script only builds the --input samplesheet. You then launch viralrecon
yourself, e.g.:

  illumina:  nextflow run nf-core/viralrecon --platform illumina \
               --input <sheet.csv> ...your other params (genome, primers, etc.)

  nanopore:  nextflow run nf-core/viralrecon --platform nanopore \
               --input <sheet.csv> --fastq_dir <dir> \
               --sequencing_summary <file> ...your other params

The trailing "..." just means the rest of your viralrecon parameters, which have
nothing to do with building the samplesheet.

Usage:
    # Illumina, a whole directory (no --platform needed; it is the default)
    ./make_viralrecon_samplesheet.py reads/ samplesheet.csv

    # Illumina from a Google Cloud Storage bucket
    ./make_viralrecon_samplesheet.py gs://my-bucket/reads/ samplesheet.csv

    # Illumina, a glob of specific files (quote it)
    ./make_viralrecon_samplesheet.py 'reads/*.fastq.gz' samplesheet.csv

    # Illumina, one sample out of a directory of many
    ./make_viralrecon_samplesheet.py reads/ one.csv --include 'SampleA*'

    # Nanopore (point at the fastq_pass directory)
    ./make_viralrecon_samplesheet.py fastq_pass/ samplesheet.csv --platform nanopore

    # Nanopore with real sample names
    ./make_viralrecon_samplesheet.py fastq_pass/ samplesheet.csv \
        --platform nanopore --barcode-map names.csv
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
FASTQ_EXTS = (".fastq.gz", ".fq.gz")
BARCODE_RE = re.compile(r"^barcode(\d+)$", re.IGNORECASE)


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
    subdirectories are not searched - or a shell glob like 'reads/*.fastq.gz'.
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


def list_barcode_dirs(location):
    """Map barcode_number(int) -> directory name for immediate barcodeNN subdirs."""
    if location.startswith("gs://"):
        # Immediate children only; gs:// subdirectories come back as prefixes ending /.
        names = [ln.strip().rstrip("/").split("/")[-1]
                 for ln in run_gsutil_ls(location.rstrip("/") + "/")
                 if ln.strip().endswith("/")]
    else:
        names = [p.name for p in Path(location).iterdir() if p.is_dir()]
    barcodes = {}
    for name in names:
        match = BARCODE_RE.match(name)
        if match:
            barcodes[int(match.group(1))] = name
    return barcodes


def pair_short_reads(files):
    """Group paired/single short reads by sample into r1/r2 slots."""
    samples = defaultdict(dict)
    for uri in files:
        name = uri.split("/")[-1]
        match = PAIR_RE.match(name)
        if not match:
            print(f"WARNING: '{name}' has no R1/R2 read number; skipping", file=sys.stderr)
            continue
        samples[match.group("sample")][f"r{match.group('read')}"] = uri
    return samples


def sanitize_sample(name):
    """viralrecon converts spaces in sample names to underscores; do it up front."""
    if " " in name:
        fixed = name.replace(" ", "_")
        print(f"WARNING: sample '{name}' contains spaces; using '{fixed}'", file=sys.stderr)
        return fixed
    return name


def load_barcode_map(path):
    """barcode_number(int) -> sample name, from a 2-column CSV/TSV 'barcode,sample'."""
    mapping = {}
    delim = "\t" if path.endswith(".tsv") else ","
    try:
        lines = Path(path).read_text().splitlines()
    except OSError as exc:
        sys.exit(f"cannot read --barcode-map {path}: {exc}")
    for n, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.lower().startswith(("barcode,", "sample,")):
            continue  # skip blank lines and an optional header row
        parts = [p.strip() for p in line.split(delim)]
        if len(parts) < 2:
            sys.exit(f"--barcode-map line {n}: expected 2 columns 'barcode,sample'")
        try:
            bc = int(re.sub(r"^barcode", "", parts[0], flags=re.IGNORECASE))
        except ValueError:
            sys.exit(f"--barcode-map line {n}: barcode '{parts[0]}' is not a number")
        mapping[bc] = parts[1]
    return mapping


def build_illumina(args):
    """Return (header, rows) for an Illumina sample,fastq_1,fastq_2 sheet."""
    files = resolve_inputs(args.input_location, args.include)
    if not files:
        sys.exit(f"No FASTQ files found at {args.input_location}")
    rows = []
    for sample, reads in sorted(pair_short_reads(files).items()):
        r1 = reads.get("r1", "")
        r2 = reads.get("r2", "")  # empty for single-end
        if not r1:
            print(f"WARNING: sample '{sample}' has R2 but no R1; skipping", file=sys.stderr)
            continue
        rows.append([sanitize_sample(sample), r1, r2])
    return ["sample", "fastq_1", "fastq_2"], rows


def build_nanopore(args):
    """Return (header, rows) for a Nanopore sample,barcode sheet."""
    barcodes = list_barcode_dirs(args.input_location)
    if not barcodes:
        sys.exit(f"No barcodeNN directories found in {args.input_location}")
    bmap = load_barcode_map(args.barcode_map) if args.barcode_map else {}
    rows = []
    for num in sorted(barcodes):
        sample = sanitize_sample(bmap.get(num, barcodes[num]))
        rows.append([sample, str(num)])
    return ["sample", "barcode"], rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # --- Required positionals (every run) ---
    parser.add_argument(
        "input_location",
        help="REQUIRED. illumina: a directory of FASTQs (non-recursive) or a quoted "
             "glob like 'reads/*.fastq.gz'. nanopore: the fastq_dir of barcodeNN/ folders.")
    parser.add_argument(
        "output", help="REQUIRED. Output samplesheet path (CSV).")

    # --- Optional ---
    parser.add_argument(
        "--platform", choices=["illumina", "nanopore"], default="illumina",
        help="OPTIONAL. 'illumina' (default) or 'nanopore'. Omit for Illumina.")
    parser.add_argument(
        "--include", action="append", metavar="GLOB",
        help="OPTIONAL (illumina only). Keep only files whose basename matches this "
             "glob, e.g. 'SampleA*'. Repeatable. Quote it.")
    parser.add_argument(
        "--barcode-map", default=None, metavar="FILE",
        help="OPTIONAL (nanopore only). 2-column CSV/TSV 'barcode,sample' giving a "
             "sample name per barcode; unlisted barcodes use their folder name.")

    args = parser.parse_args()

    # Reject options paired with the wrong platform rather than silently ignoring them.
    if args.platform == "illumina" and args.barcode_map:
        sys.exit("--barcode-map only applies to --platform nanopore")
    if args.platform == "nanopore" and args.include:
        sys.exit("--include only applies to --platform illumina")

    if args.platform == "illumina":
        header, rows = build_illumina(args)
    else:
        header, rows = build_nanopore(args)

    if not rows:
        sys.exit("No samples assembled; check your inputs")

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")  # comma: viralrecon requires CSV
        writer.writerow(header)
        writer.writerows(rows)

    print(f"Wrote {len(rows)} samples to {args.output} (--platform {args.platform})")


if __name__ == "__main__":
    main()
