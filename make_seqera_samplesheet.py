#!/usr/bin/env python3

import argparse
import csv
import fnmatch
import sys
from collections import defaultdict
from pathlib import Path
import re


PAIR_RE = re.compile(r"^(?P<sample>.+?)_R?(?P<read>[12])(?:_\d{3})?\.(?:fastq|fq)\.gz$")
FASTQ_EXTS = (".fastq.gz", ".fq.gz")


def read_filenames(metadata_csv):
    """Return the `filename` values from the metadata CSV, in file order.

    Errors clearly if the file is missing or has no `filename` column. Blank
    filename cells are skipped.
    """
    try:
        text = Path(metadata_csv).read_text()
    except OSError as exc:
        sys.exit(f"cannot read metadata CSV {metadata_csv}: {exc}")
    reader = csv.DictReader(text.splitlines())
    if reader.fieldnames is None:
        sys.exit(f"metadata CSV {metadata_csv} is empty")
        
    field = next((f for f in reader.fieldnames if f.strip().lower() == "filename"), None)
    if field is None:
        sys.exit(f"metadata CSV {metadata_csv} has no 'filename' column "
                 f"(header: {', '.join(reader.fieldnames)})")
    return [row[field].strip() for row in reader if (row.get(field) or "").strip()]


def select_fastqs(filenames, patterns=None):
    """Keep only FASTQ filenames, warning on the rest; apply optional --include."""
    fastqs = []
    for name in filenames:
        if not name.endswith(FASTQ_EXTS):
            print(f"WARNING: '{name}' is not a FASTQ (.fastq.gz/.fq.gz); skipping",
                  file=sys.stderr)
            continue
        if patterns and not any(fnmatch.fnmatch(name, pat) for pat in patterns):
            continue
        fastqs.append(name)
    return fastqs


def pair_short_reads(filenames):
    """Group paired/single short reads by sample into r1/r2 slots."""
    samples = defaultdict(dict)
    for name in filenames:
        match = PAIR_RE.match(name)
        if not match:
            print(f"WARNING: '{name}' has no R1/R2 read number; skipping", file=sys.stderr)
            continue
        samples[match.group("sample")][f"r{match.group('read')}"] = name
    return samples


def sanitize_sample(name):
    """viralrecon converts spaces in sample names to underscores; do it up front."""
    if " " in name:
        fixed = name.replace(" ", "_")
        print(f"WARNING: sample '{name}' contains spaces; using '{fixed}'", file=sys.stderr)
        return fixed
    return name


def with_prefix(prefix, filename):
    """Join --prefix and a filename with exactly one '/'. Empty prefix -> bare name."""
    if not prefix:
        return filename
    return prefix.rstrip("/") + "/" + filename


def build_rows(filenames, prefix):
    """Return (header, rows) for a sample,fastq_1,fastq_2 sheet."""
    rows = []
    for sample, reads in sorted(pair_short_reads(filenames).items()):
        r1 = reads.get("r1")
        r2 = reads.get("r2")  # may be None for single-end
        if not r1:
            print(f"WARNING: sample '{sample}' has R2 but no R1; skipping", file=sys.stderr)
            continue
        rows.append([
            sanitize_sample(sample),
            with_prefix(prefix, r1),
            with_prefix(prefix, r2) if r2 else "",
        ])
    return ["sample", "fastq_1", "fastq_2"], rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "metadata_csv",
        help="REQUIRED. The APGAP metadata CSV (must have a 'filename' column).")
    parser.add_argument(
        "output", help="REQUIRED. Output samplesheet path (CSV).")
    parser.add_argument(
        "--prefix", default="", metavar="PREFIX",
        help="OPTIONAL but recommended. Bucket/dir prepended to each filename to form "
             "the full path, e.g. 'gs://my-bucket/reads/'. Omit to write bare filenames.")
    parser.add_argument(
        "--include", action="append", metavar="GLOB",
        help="OPTIONAL. Keep only filenames whose basename matches this glob, e.g. "
             "'SampleA*'. Repeatable. Quote it.")

    args = parser.parse_args()

    if not args.prefix:
        print("WARNING: no --prefix given; writing bare filenames. The sheet will not "
              "run on Google Cloud Batch without resolvable paths.", file=sys.stderr)

    filenames = read_filenames(args.metadata_csv)
    fastqs = select_fastqs(filenames, args.include)
    if not fastqs:
        sys.exit(f"No FASTQ filenames found in {args.metadata_csv}")

    header, rows = build_rows(fastqs, args.prefix)
    if not rows:
        sys.exit("No samples assembled; check your inputs")

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")  # comma: nf-core requires CSV
        writer.writerow(header)
        writer.writerows(rows)

    print(f"Wrote {len(rows)} samples to {args.output}")


if __name__ == "__main__":
    main()
