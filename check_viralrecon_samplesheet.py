#!/usr/bin/env python3
"""
Validate an nf-core/viralrecon samplesheet (Illumina or Nanopore).

Auto-detects the platform from the header and applies the matching rules:
  Illumina : sample,fastq_1,fastq_2   (fastq_2 empty = single-end)
  Nanopore : sample,barcode           (barcode must be a positive integer)

viralrecon requires comma-separated input. Exits 1 on any ERROR, 0 if clean
(WARNINGs do not fail).

Checks
  Common (ERROR):
    - file missing / empty / contains carriage returns (CRLF / ^M)
    - file looks tab-delimited (viralrecon requires CSV)
    - header is not 'sample,fastq_1,fastq_2' or 'sample,barcode'
    - --platform given but the header is for the other platform
    - a row whose column count does not match the header
    - empty sample name
  Common (WARNING):
    - missing trailing newline
    - sample name with spaces (viralrecon converts them to underscores) or with
      characters outside letters/numbers/._-

  Illumina (ERROR):
    - fastq_1 empty, or fastq_1/fastq_2 not ending in .fastq.gz / .fq.gz
    - a fully identical duplicate row
  Illumina (WARNING):
    - the same FASTQ path reused on more than one row
    - repeated sample names (allowed; viralrecon merges them as replicates)

  Nanopore (ERROR):
    - barcode not a positive integer
    - duplicate barcode
  Nanopore (WARNING):
    - the same sample name on more than one barcode

NOTE: duplicate sample names are NOT an error for Illumina - viralrecon merges
same-named rows as technical replicates. Read-file existence is not checked.

Usage:
    ./check_viralrecon_samplesheet.py samplesheet.csv
    ./check_viralrecon_samplesheet.py samplesheet.csv --platform nanopore
"""
import argparse
import re
import sys

ILLUMINA_HEADER = ["sample", "fastq_1", "fastq_2"]
NANOPORE_HEADER = ["sample", "barcode"]
FASTQ_EXTS = (".fastq.gz", ".fq.gz")
SAMPLE_SAFE = re.compile(r"^[A-Za-z0-9._-]+$")
INT_RE = re.compile(r"^\d+$")


def detect_platform(header_line):
    """Return (platform, note). platform is 'illumina'/'nanopore' or None.
    When None, note is 'tab' (right columns but tab-delimited) or 'unknown'."""
    comma = [h.strip() for h in header_line.split(",")]
    if comma == ILLUMINA_HEADER:
        return "illumina", ""
    if comma == NANOPORE_HEADER:
        return "nanopore", ""
    tab = [h.strip() for h in header_line.split("\t")]
    if tab in (ILLUMINA_HEADER, NANOPORE_HEADER):
        return None, "tab"
    return None, "unknown"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("samplesheet", help="Path to the viralrecon samplesheet (CSV)")
    parser.add_argument("--platform", choices=["illumina", "nanopore"], default=None,
                        help="OPTIONAL. Assert the sheet is for this platform and error on "
                             "mismatch. Default: auto-detect from the header.")
    args = parser.parse_args()

    errors = []
    warnings = []

    # --- File-level checks ---
    try:
        with open(args.samplesheet, "rb") as fh:
            raw = fh.read()
    except OSError as exc:
        sys.exit(f"ERROR: cannot read {args.samplesheet}: {exc}")

    if not raw.strip():
        sys.exit(f"ERROR: {args.samplesheet} is empty")

    if b"\r" in raw:
        bad = [i + 1 for i, ln in enumerate(raw.split(b"\n")) if ln.endswith(b"\r")]
        preview = ", ".join(map(str, bad[:10])) + (" ..." if len(bad) > 10 else "")
        errors.append(f"carriage returns (CRLF / ^M) on line(s): {preview}. Fix with: "
                      f"tr -d '\\r' < {args.samplesheet} > fixed && mv fixed {args.samplesheet}")

    if not raw.endswith(b"\n"):
        warnings.append("file has no trailing newline")

    lines = raw.decode("utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        sys.exit(f"ERROR: {args.samplesheet} has a header but no sample rows")

    # --- Platform / header ---
    platform, note = detect_platform(lines[0])
    if platform is None and note == "tab":
        sys.exit(f"ERROR: {args.samplesheet} looks tab-delimited; viralrecon requires "
                 "comma-separated (CSV). Re-save with commas.")
    if platform is None:
        sys.exit(f"ERROR: unrecognized header {[h.strip() for h in lines[0].split(',')]}; "
                 f"expected {ILLUMINA_HEADER} or {NANOPORE_HEADER}")
    if args.platform and args.platform != platform:
        sys.exit(f"ERROR: --platform {args.platform} requested but the header is a "
                 f"{platform} sheet")

    ncol = len(ILLUMINA_HEADER if platform == "illumina" else NANOPORE_HEADER)

    seen_samples = {}    # sample -> first line
    seen_barcodes = {}   # nanopore: barcode -> first line
    seen_rows = set()    # illumina: full-row dedupe
    seen_fastqs = {}     # illumina: fastq path -> first line
    replicate_samples = set()
    data_rows = 0

    for n, line in enumerate(lines[1:], start=2):
        if line.strip() == "":
            warnings.append(f"line {n}: blank line skipped")
            continue
        fields = [f.strip() for f in line.split(",")]
        if len(fields) != ncol:
            errors.append(f"line {n}: expected {ncol} columns, found {len(fields)}")
            continue
        data_rows += 1

        # Common: sample name
        sample = fields[0]
        if not sample:
            errors.append(f"line {n}: empty sample name")
        elif " " in sample:
            warnings.append(f"line {n}: sample '{sample}' has spaces "
                            "(viralrecon will convert them to underscores)")
        elif not SAMPLE_SAFE.match(sample):
            warnings.append(f"line {n}: sample '{sample}' has unusual characters; "
                            "stick to letters/numbers/._-")

        if platform == "illumina":
            if tuple(fields) in seen_rows:
                errors.append(f"line {n}: duplicate row")
                continue
            seen_rows.add(tuple(fields))

            fq1, fq2 = fields[1], fields[2]
            if not fq1:
                errors.append(f"line {n}: fastq_1 is empty (required)")
            elif not fq1.endswith(FASTQ_EXTS):
                errors.append(f"line {n}: fastq_1 '{fq1}' must end in .fastq.gz/.fq.gz")
            if fq2 and not fq2.endswith(FASTQ_EXTS):
                errors.append(f"line {n}: fastq_2 '{fq2}' must end in .fastq.gz/.fq.gz "
                              "(leave it empty for single-end)")

            for fq in (fq1, fq2):
                if fq:
                    if fq in seen_fastqs:
                        warnings.append(f"line {n}: FASTQ '{fq}' also used on line {seen_fastqs[fq]}")
                    else:
                        seen_fastqs[fq] = n

            if sample:
                if sample in seen_samples:
                    replicate_samples.add(sample)
                else:
                    seen_samples[sample] = n

        else:  # nanopore
            barcode = fields[1]
            if not INT_RE.match(barcode) or int(barcode) < 1:
                errors.append(f"line {n}: barcode '{barcode}' must be a positive integer")
            elif barcode in seen_barcodes:
                errors.append(f"line {n}: duplicate barcode '{barcode}' "
                              f"(also line {seen_barcodes[barcode]})")
            else:
                seen_barcodes[barcode] = n

            if sample:
                if sample in seen_samples:
                    warnings.append(f"line {n}: sample '{sample}' repeats "
                                    f"(also line {seen_samples[sample]}); expected one per barcode")
                else:
                    seen_samples[sample] = n

    if replicate_samples:
        warnings.append(f"repeated sample name(s) {sorted(replicate_samples)} will be merged "
                        "as technical replicates (allowed)")

    # --- Report ---
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)

    if errors:
        print(f"FAILED: {len(errors)} error(s), {len(warnings)} warning(s) in {args.samplesheet}",
              file=sys.stderr)
        sys.exit(1)

    print(f"OK: {platform} samplesheet, {data_rows} row(s) valid, "
          f"{len(warnings)} warning(s) in {args.samplesheet}")


if __name__ == "__main__":
    main()
