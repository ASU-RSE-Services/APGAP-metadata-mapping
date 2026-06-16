#!/usr/bin/env python3
"""
Validate a Bactopia FOFN (samplesheet). Works for v3.x and v4.x - the 7-column
format is shared.

Format (TAB-delimited, 7 columns):
    sample  runtype  genome_size  species  r1  r2  extra

Column rules depend on each row's runtype, and a single FOFN may MIX runtypes
(Bactopia processes every sample by its own runtype - that is allowed, not an
error). Exits 1 on any ERROR, 0 if clean (WARNINGs do not fail).

  runtype       r1                      r2            extra
  paired-end    R1 (.fastq.gz/.fq.gz)   R2 fastq      empty
  single-end    read fastq              empty         empty
  ont           ONT read fastq          empty         empty
  hybrid        R1 fastq                R2 fastq      long read fastq
  short_polish  R1 fastq                R2 fastq      long read fastq
  assembly      FASTA (.fna/.fasta/.fa .gz)  empty    empty

Checks
  Common (ERROR):
    - file missing / empty / contains carriage returns (CRLF / ^M)
    - file looks comma-delimited (the FOFN must be tab-delimited)
    - header is not 'sample runtype genome_size species r1 r2 extra'
    - a row without exactly 7 columns
    - empty sample name, or a duplicate sample name (Bactopia requires unique)
    - runtype not one of the six valid values
    - genome_size present but not a non-negative integer
    - read/assembly columns that violate the row's runtype rules above
  Common (WARNING):
    - missing trailing newline
    - sample name with spaces or characters outside letters/numbers/._-
    - empty genome_size (Bactopia prepare writes 0) or empty species
      (UNKNOWN_SPECIES is the default)

File existence of the reads is not checked. Usage:
    ./check_bactopia_samplesheet.py samples.tsv
"""
import argparse
import re
import sys

FOFN_HEADER = ["sample", "runtype", "genome_size", "species", "r1", "r2", "extra"]
FASTQ_EXTS = (".fastq.gz", ".fq.gz")
FASTA_EXTS = (".fna.gz", ".fasta.gz", ".fa.gz")
VALID_RUNTYPES = {"paired-end", "single-end", "ont", "hybrid", "short_polish", "assembly"}
SAMPLE_SAFE = re.compile(r"^[A-Za-z0-9._-]+$")
INT_RE = re.compile(r"^\d+$")


def check_read_columns(runtype, r1, r2, extra, n, errors):
    """Apply the per-runtype rules for the r1/r2/extra columns."""
    def present(v):
        return v != ""

    def fastq(v):
        return v.endswith(FASTQ_EXTS)

    if runtype in ("paired-end", "hybrid", "short_polish"):
        if not present(r1) or not fastq(r1):
            errors.append(f"line {n}: {runtype} needs r1 as a .fastq.gz/.fq.gz file")
        if not present(r2) or not fastq(r2):
            errors.append(f"line {n}: {runtype} needs r2 as a .fastq.gz/.fq.gz file")
        if runtype == "paired-end":
            if present(extra):
                errors.append(f"line {n}: paired-end must leave extra empty")
        elif not present(extra) or not fastq(extra):
            errors.append(f"line {n}: {runtype} needs extra as the long-read "
                          ".fastq.gz/.fq.gz file")

    elif runtype in ("single-end", "ont"):
        if not present(r1) or not fastq(r1):
            errors.append(f"line {n}: {runtype} needs r1 as a .fastq.gz/.fq.gz file")
        if present(r2):
            errors.append(f"line {n}: {runtype} must leave r2 empty")
        if present(extra):
            errors.append(f"line {n}: {runtype} must leave extra empty")

    elif runtype == "assembly":
        if not present(r1) or not r1.endswith(FASTA_EXTS):
            errors.append(f"line {n}: assembly needs r1 as a FASTA "
                          "(.fna.gz/.fasta.gz/.fa.gz)")
        if present(r2):
            errors.append(f"line {n}: assembly must leave r2 empty")
        if present(extra):
            errors.append(f"line {n}: assembly must leave extra empty")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("samplesheet", help="Path to the Bactopia FOFN (tab-delimited)")
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

    # --- Delimiter / header ---
    header = [h.strip() for h in lines[0].split("\t")]
    if header != FOFN_HEADER:
        if [h.strip() for h in lines[0].split(",")] == FOFN_HEADER:
            sys.exit(f"ERROR: {args.samplesheet} looks comma-delimited; the Bactopia "
                     "FOFN must be tab-delimited. Re-save with tabs.")
        sys.exit(f"ERROR: header must be {FOFN_HEADER}, got {header}")

    # --- Per-row checks ---
    seen = {}
    for n, line in enumerate(lines[1:], start=2):
        if line.strip() == "":
            warnings.append(f"line {n}: blank line skipped")
            continue
        fields = line.split("\t")
        if len(fields) != len(FOFN_HEADER):
            errors.append(f"line {n}: expected {len(FOFN_HEADER)} tab-separated columns, "
                          f"found {len(fields)}")
            continue
        sample, runtype, gsize, species, r1, r2, extra = [f.strip() for f in fields]

        # sample
        if not sample:
            errors.append(f"line {n}: empty sample name")
        elif sample in seen:
            errors.append(f"line {n}: duplicate sample '{sample}' (also line {seen[sample]}); "
                          "Bactopia requires unique sample names")
        else:
            seen[sample] = n
            if " " in sample:
                warnings.append(f"line {n}: sample '{sample}' has spaces")
            elif not SAMPLE_SAFE.match(sample):
                warnings.append(f"line {n}: sample '{sample}' has unusual characters; "
                                "stick to letters/numbers/._-")

        # runtype
        if runtype not in VALID_RUNTYPES:
            errors.append(f"line {n}: runtype '{runtype}' is not one of "
                          f"{sorted(VALID_RUNTYPES)}")
        else:
            check_read_columns(runtype, r1, r2, extra, n, errors)

        # genome_size
        if gsize == "":
            warnings.append(f"line {n}: empty genome_size (Bactopia prepare writes 0)")
        elif not INT_RE.match(gsize):
            errors.append(f"line {n}: genome_size '{gsize}' must be a non-negative integer")

        # species
        if species == "":
            warnings.append(f"line {n}: empty species (UNKNOWN_SPECIES is the default)")

    # --- Report ---
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)

    if errors:
        print(f"FAILED: {len(errors)} error(s), {len(warnings)} warning(s) in {args.samplesheet}",
              file=sys.stderr)
        sys.exit(1)

    print(f"OK: {len(seen)} valid sample(s), {len(warnings)} warning(s) in {args.samplesheet}")


if __name__ == "__main__":
    main()
