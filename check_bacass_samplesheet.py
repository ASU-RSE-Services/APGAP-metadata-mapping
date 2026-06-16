#!/usr/bin/env python3
"""
Validate an nf-core/bacass samplesheet.

Catches structural and field problems here instead of after a pipeline launch.
Exits 1 if any ERROR is found, 0 if the sheet is clean (WARNINGs do not fail).

Expected format (6 columns): ID, R1, R2, LongFastQ, Fast5, GenomeSize
  Delimiter follows the file extension, exactly as bacass reads it:
  .csv -> comma, .tsv -> tab.

Checks performed
  ERROR (fail the check):
    - file missing, empty, or containing carriage returns (CRLF / ^M)
    - content delimiter disagrees with the extension (bacass would mis-parse)
    - header is not exactly ID,R1,R2,LongFastQ,Fast5,GenomeSize
    - a row does not have exactly 6 columns
    - ID empty, not starting with a letter, or with characters other than
      letters / numbers / underscores
    - duplicate ID
    - an R1/R2/LongFastQ value that is empty, or non-NA but not a .fastq.gz /
      .fq.gz file
    - R2 present without R1
    - a row with no reads at all (R1, R2, LongFastQ all NA)
    - GenomeSize that is non-NA but not a number ending in 'm' (e.g. 4.6m)
  WARNING (does not fail):
    - missing trailing newline
    - header names correct but not in canonical order
    - R1 present but R2=NA (single-end short reads are unusual for bacass)
    - GenomeSize=NA on a long-read or hybrid row (assemblers want the estimate)
    - rows mix assembly modes (short / long / hybrid); one --assembly_type
      cannot fit all rows

Does NOT check whether the read files actually exist (local and gs:// paths are
accepted as written).

Usage:
    ./check_bacass_samplesheet.py samplesheet.csv
"""
import argparse
import re
import sys

EXPECTED_HEADER = ["ID", "R1", "R2", "LongFastQ", "Fast5", "GenomeSize"]
FASTQ_EXTS = (".fastq.gz", ".fq.gz")
ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
GENOME_RE = re.compile(r"^\d+(?:\.\d+)?m$")


def dname(delim):
    return "tab" if delim == "\t" else "comma"


def present(value):
    """A read column counts as present only if it is not empty and not NA."""
    return value not in ("", "NA")


def detect_delimiter(header_line, path):
    """Return (delimiter bacass infers from the extension, delimiter the content implies)."""
    by_ext = "\t" if path.endswith(".tsv") else ","
    comma_n = len(header_line.split(","))
    tab_n = len(header_line.split("\t"))
    want = len(EXPECTED_HEADER)
    if tab_n == want and comma_n != want:
        by_content = "\t"
    elif comma_n == want and tab_n != want:
        by_content = ","
    else:
        by_content = by_ext  # ambiguous; trust the extension
    return by_ext, by_content


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("samplesheet", help="Path to the bacass samplesheet (.csv or .tsv)")
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
        errors.append(
            f"carriage returns (CRLF / ^M) found on line(s): {preview}. Fix with: "
            f"tr -d '\\r' < {args.samplesheet} > fixed && mv fixed {args.samplesheet}")

    if not raw.endswith(b"\n"):
        warnings.append("file has no trailing newline")

    lines = raw.decode("utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        sys.exit(f"ERROR: {args.samplesheet} has a header but no sample rows")

    # --- Delimiter / header ---
    by_ext, by_content = detect_delimiter(lines[0], args.samplesheet)
    if by_ext != by_content:
        ext = ".tsv" if args.samplesheet.endswith(".tsv") else ".csv"
        errors.append(
            f"content looks {dname(by_content)}-delimited but the {ext} extension makes "
            f"bacass read it as {dname(by_ext)}-delimited; rename the file or fix the delimiter")
    delim = by_content

    header = [h.strip() for h in lines[0].split(delim)]
    if header == EXPECTED_HEADER:
        idx = {name: i for i, name in enumerate(EXPECTED_HEADER)}
    elif sorted(header) == sorted(EXPECTED_HEADER):
        warnings.append(f"header names correct but not in canonical order: {header}")
        idx = {name: header.index(name) for name in EXPECTED_HEADER}
    else:
        errors.append(f"header must be {EXPECTED_HEADER}, got {header}")
        idx = {name: i for i, name in enumerate(EXPECTED_HEADER)}  # best-effort positional

    # --- Per-row checks ---
    seen_ids = {}
    modes = set()
    for n, line in enumerate(lines[1:], start=2):  # line 1 is the header
        if line.strip() == "":
            warnings.append(f"line {n}: blank line skipped")
            continue

        fields = [f.strip() for f in line.split(delim)]
        if len(fields) != len(EXPECTED_HEADER):
            errors.append(f"line {n}: expected {len(EXPECTED_HEADER)} columns, found {len(fields)}")
            continue

        sid = fields[idx["ID"]]
        r1 = fields[idx["R1"]]
        r2 = fields[idx["R2"]]
        lng = fields[idx["LongFastQ"]]
        gsize = fields[idx["GenomeSize"]]

        # ID
        if not sid:
            errors.append(f"line {n}: empty ID")
        elif not ID_RE.match(sid):
            errors.append(f"line {n}: ID '{sid}' must start with a letter and contain only "
                          "letters/numbers/underscores")
        elif sid in seen_ids:
            errors.append(f"line {n}: duplicate ID '{sid}' (also on line {seen_ids[sid]})")
        else:
            seen_ids[sid] = n

        # Read columns: empty -> error, else NA or a FASTQ file
        for col, val in (("R1", r1), ("R2", r2), ("LongFastQ", lng)):
            if val == "":
                errors.append(f"line {n}: {col} is empty; use NA for an absent read")
            elif val != "NA" and not val.endswith(FASTQ_EXTS):
                errors.append(f"line {n}: {col} '{val}' is not NA and not a .fastq.gz/.fq.gz file")

        has_r1, has_r2, has_long = present(r1), present(r2), present(lng)

        if has_r2 and not has_r1:
            errors.append(f"line {n}: R2 present but R1=NA (R1 is required for paired-end reads)")
        if not has_r1 and not has_r2 and not has_long:
            errors.append(f"line {n}: no reads - R1, R2 and LongFastQ are all NA")
        if has_r1 and not has_r2 and not has_long:
            warnings.append(f"line {n}: R1 present but R2=NA (single-end short reads are unusual "
                            "for bacass)")

        # GenomeSize
        if gsize != "NA" and not GENOME_RE.match(gsize):
            errors.append(f"line {n}: GenomeSize '{gsize}' must be NA or a number ending in 'm' "
                          "(e.g. 4.6m)")

        # Assembly mode for this row (for cross-row consistency)
        if has_r1 and has_long:
            mode = "hybrid"
        elif has_long and not has_r1:
            mode = "long"
        elif has_r1 and not has_long:
            mode = "short"
        else:
            mode = "unknown"
        modes.add(mode)

        if mode in ("long", "hybrid") and gsize == "NA":
            warnings.append(f"line {n}: GenomeSize=NA on a {mode} row; long-read assembly "
                            "usually wants a size estimate")

    real_modes = modes - {"unknown"}
    if len(real_modes) > 1:
        warnings.append(f"sheet mixes assembly modes {sorted(real_modes)}; a single "
                        "--assembly_type cannot fit all rows (bacass runs every sample with one)")

    # --- Report ---
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)

    if errors:
        print(f"FAILED: {len(errors)} error(s), {len(warnings)} warning(s) "
              f"in {args.samplesheet}", file=sys.stderr)
        sys.exit(1)

    inferred = next(iter(real_modes)) if len(real_modes) == 1 else "mixed"
    suffix = f"; --assembly_type {inferred}" if inferred != "mixed" else ""
    print(f"OK: {len(seen_ids)} valid sample(s), {len(warnings)} warning(s) "
          f"in {args.samplesheet}{suffix}")


if __name__ == "__main__":
    main()
