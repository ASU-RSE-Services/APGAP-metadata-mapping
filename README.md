# APGAP-metadata-mapping

Repository for scripts that map from APGAP metadata to other metadata formats.

## Scripts

### `make_seqera_samplesheet.py` — metadata CSV → Seqera samplesheet

Converts an **APGAP metadata CSV** (the `_apgap_metadata.csv` export, whose first
column is `filename`) into a Seqera/Nextflow Illumina samplesheet with columns
`sample,fastq_1,fastq_2`.

It reads the `filename` column (other columns are ignored), keeps only FASTQ rows
(`.fastq.gz` / `.fq.gz`), pairs R1/R2 by filename, and prepends a bucket `--prefix`
to each filename to form the full path. A sample with only R1 becomes single-end
(empty `fastq_2`).

```bash
# reads live in a gs:// bucket
./make_seqera_samplesheet.py _apgap_metadata.csv samplesheet.csv \
    --prefix gs://my-bucket/reads/

# just one sample out of the csv
./make_seqera_samplesheet.py _apgap_metadata.csv one.csv \
    --prefix gs://my-bucket/reads/ --include 'SampleA*'
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `metadata_csv` | yes | the APGAP metadata CSV (must have a `filename` column) |
| `output` | yes | output samplesheet path (written as CSV) |
| `--prefix` | recommended | bucket/dir prepended to each filename, e.g. `gs://my-bucket/reads/`; omit to write bare filenames (warns) |
| `--include` | no | keep only filenames matching this glob, e.g. `'SampleA*'`; repeatable |

The output format matches `make_viralrecon_samplesheet.py`'s Illumina sheet, so it
can be validated with `check_viralrecon_samplesheet.py samplesheet.csv`.

### Pipeline samplesheet builders (scan a FASTQ directory)

These build a samplesheet by scanning a directory/glob of FASTQ files (local or
`gs://`), rather than from a metadata CSV:

- `make_viralrecon_samplesheet.py` — nf-core/viralrecon (Illumina or Nanopore).
- `make_bacass_samplesheet.py` — nf-core/bacass (short / long / hybrid).
- `make_bactopia_fofn.py` — Bactopia FOFN (v3.x / v4.x).

### Validators

- `check_viralrecon_samplesheet.py`
- `check_bacass_samplesheet.py`
- `check_bactopia_samplesheet.py`

Run any script with `--help` for full usage.
