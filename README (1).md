# GWAMA SNP liftover

A reproducible GRCh37-to-GRCh38 conversion workflow for **19-column GWAMA-style summary statistics** whose chromosome, position and allele labels are embedded in `rs_number`, for example `10:100000625_A_G`.

The workflow uses **UCSC liftOver as the coordinate engine**. Its contribution is the surrounding format-specific parsing, effect-allele bookkeeping, GRCh38 SNP validation, review records, streaming execution, resumable checkpoints and final integrity audit. It is not a new coordinate-mapping algorithm.

The original application processed European-ancestry MCH, MCV and RDW summary statistics. The portable scripts also accept other named traits **with the same input schema and GRCh37 coordinate convention**.

## Scope

- Accepted output contains uniquely mapped, canonical-chromosome, biallelic **SNPs** whose transformed allele pair includes the GRCh38 FASTA base.
- Every other input row is retained in a review file with a reason. Non-SNP or unsupported-allele records require separate variant-aware conversion; they are not certified GRCh38 variants here.
- Statistical values are preserved. This stage does not apply new genomic control, choose a primary P-value, filter on allele frequency or sample size, perform LD harmonisation, globally deduplicate variants, fine-map loci or integrate a GRN.
- Full GWAS files, FASTA references, UCSC executables and chain files are not distributed in this repository.

## Why the allele fields matter

In GWAMA output, **`reference_allele` is the effect allele and `other_allele` is the non-effect allele**, as documented by [GWAMA](https://genomics.ut.ee/en/tools). The effect allele need not be the assembly reference base. The two alleles within the composite marker are used to check the allele pair, without assuming that the first label is the genome REF.

This workflow maintains two separate representations:

| Fields | Meaning |
|---|---|
| `reference_allele_GRCh37`, `other_allele_GRCh37` | Original reported effect and non-effect alleles |
| `reference_allele_hg38`, `other_allele_hg38` | The same allele identities expressed on the target forward strand |
| `REF_hg38`, `ALT_hg38` | Assembly reference base and the other allele, assigned using GRCh38 FASTA |
| `ALT_hg38_from` | Which original reported allele corresponds to genomic ALT |

For reverse-strand mappings both reported SNP alleles are complemented. Their identities are not swapped, so **`beta`, `z`, `eaf`, confidence limits, P-values and all other statistical values remain unchanged**. The reported effect in `beta` refers to `reference_allele_hg38`; it is **not automatically an ALT-relative effect**. If an application later requires effects relative to genomic ALT, it must explicitly use `ALT_hg38_from` to determine any sign/frequency transformation.

## Requirements

- Python **3.9 or later** on Linux; the workflow uses the Python standard library, including POSIX file locking and memory mapping.
- A UCSC **liftOver** executable compatible with the host architecture.
- The official **`hg19ToHg38.over.chain.gz`** forward chain.
- An **uncompressed GRCh38 primary-assembly FASTA** and its matching `.fai`, with all of chromosomes 1–22, X and Y. UCSC-style names and unprefixed canonical names are supported in the FASTA index.
- `pigz` is optional; Python gzip is the fallback. `samtools` is needed only if the FASTA index must be generated.

No pandas, NumPy, GWASLab, bcftools or R packages are required for this SNP workflow.

## Input schema

The original files are named:

```text
MCH_EA_GWAMA_20190214_GC.gz
MCV_EA_GWAMA_20190214_GC.gz
RDW_EA_GWAMA_20190214_GC.gz
```

They are gzip-compressed text files with the following 19 columns:

```text
rs_number reference_allele other_allele eaf beta se beta_95L beta_95U z p-value _-log10_p-value q_statistic q_p-value i2 n_studies n_samples effects chi-squared_GCcorr p-value_GCcorr
```

Tabs are preferred; whitespace-delimited input is also handled. Column access uses header names. The production parser requires 19 uniquely named columns and the required marker, allele and statistical fields. `rs_number` must have the form `CHR:POS_A_B` or `chrCHR:POS_A_B`, with **1-based GRCh37 POS**. Human canonical chromosomes 1–22, X and Y are the supported analysis scope. The source build must be established from the dataset provenance; a composite marker alone does not establish the build. Other schemas must be explicitly adapted before using this workflow.

## Quick start

Run the following commands from the extracted repository directory.

### 1. Prepare local data and references

```bash
mkdir -p data/GWAS/raw tools reference/liftover reference/GRCh38
python3 scripts/fetch_ucsc_resources.py --binary tools/liftOver
```

The downloader retrieves the official forward chain and, when `--binary` is supplied, the **Linux x86_64** UCSC executable. Existing files are preserved. On another architecture, supply a compatible executable yourself. External-resource terms are linked in [docs/references.md](docs/references.md).

Place the three GWAS inputs under `data/GWAS/raw/` and provide a validated GRCh38 FASTA at the path configured below. If its index is absent:

```bash
samtools faidx reference/GRCh38/GRCh38.primary_assembly.genome.fa
```

The source filenames are listed by the [NHLBI GRASP full-results resource](https://www.nhlbi.nih.gov/grasp/full-results); the repository does not mirror the data. Dataset terms and original study citations remain applicable.

### 2. Configure paths

```bash
cp config.example.json config.local.json
```

Edit `config.local.json`. Paths inside a JSON profile are resolved relative to **that profile's directory**. Explicit command-line paths are resolved relative to the current working directory. Absolute paths are also accepted.

The input mapping can contain one or more safe trait labels. To override it on the command line, repeat `--input TRAIT=FILE`; these overrides replace the entire profile input mapping.

```bash
python3 scripts/Step_01_inspect.py --config config.local.json
python3 scripts/Step_02_pilot.py --config config.local.json --rows 10000
```

The pilot writes to a fresh `pilot_*` directory next to the configured production directory. It samples the **first** N rows, not random genome-wide rows. A pilot lacking reverse-strand examples does not establish reverse-strand performance.

### 3. Run production

```bash
python3 -u scripts/Step_03_liftover.py run --config config.local.json
```

Alternatively, provide all paths directly:

```bash
python3 -u scripts/Step_03_liftover.py run \
  --input MCH=data/GWAS/raw/MCH_EA_GWAMA_20190214_GC.gz \
  --input MCV=data/GWAS/raw/MCV_EA_GWAMA_20190214_GC.gz \
  --input RDW=data/GWAS/raw/RDW_EA_GWAMA_20190214_GC.gz \
  --liftover tools/liftOver \
  --chain reference/liftover/hg19ToHg38.over.chain.gz \
  --fasta reference/GRCh38/GRCh38.primary_assembly.genome.fa \
  --out results/run
```

Defaults are 500,000 rows per chunk and two pigz compression threads. Traits run sequentially. The production engine checks the forward chain direction, target FASTA chromosome lengths, input/resource fingerprints and output-directory lock before processing.

### 4. Monitor and audit

```bash
python3 scripts/Step_03_liftover.py status --config config.local.json
python3 -u scripts/Step_04_final_audit.py --config config.local.json
```

Completion markers are `ALL_TRAITS_COMPLETE` for production and **`SNP_LIFTOVER_FINAL_AUDIT_PASS`** for the final audit. The audit independently reads the final compressed files to compare SHA256 values with production records, verifies metadata fingerprints and checkpoint accounting, and checks output headers and the first accepted record. It does **not** independently recount every decompressed output row or repeat all variant mappings. Source gzip CRC checks are performed during the production scan to EOF.

## Mapping and acceptance rules

For a source SNP at 1-based position p, a BED6 interval is written as `[p−1, p)` with a persistent data-row ID, score 0 and `+` strand. The coordinate call is:

```bash
liftOver -minMatch=1.0 -multiple input.bed hg19ToHg38.over.chain.gz mapped.bed unmapped.bed
```

Repeated identical destination tuples are collapsed; multiple **distinct** chromosome/start/end/strand destinations are reviewed. A unique result must span one base on a canonical target chromosome. The target position is the BED start plus one. The transformed allele pair is compared against the GRCh38 FASTA base. An N base or an allele mismatch is reviewed, without changing a position by one base, guessing a strand or silently replacing an allele. Palindromic SNPs can pass because strand transformation comes from the chain; downstream study/LD-panel harmonisation remains separate.

The detailed schema and mutually exclusive review statuses are in [docs/input_output_schema.md](docs/input_output_schema.md).

## Outputs and restart behaviour

For each configured trait, production writes:

| File | Content |
|---|---|
| `<TRAIT>/<TRAIT>.SNP.GRCh38.tsv.gz` | Accepted SNPs, original statistical text and explicit source/target allele fields |
| `<TRAIT>/<TRAIT>.review.tsv.gz` | All other rows and their review statuses; malformed raw rows retained as base64 |
| `<TRAIT>/chunks/chunk_*/qc.json` | Chunk counters, row ranges, gzip offsets, checksums and configuration ID |
| `<TRAIT>/COMPLETE.json` | Trait-level completion, file sizes and SHA256 values |
| `QC_summary.tsv` | Trait-level counts and denominator-specific match rates |
| `run_manifest.json`, `progress.json`, `COMPLETE.json` | Provenance, current progress and overall completion |
| `final_audit/Step_04_*/` | Final audit TSV/JSON reports |

Final gzip files consist of a single header member followed by header-free compressed chunk members. They are ordinary **multi-member gzip**, not BGZF. Use a reader that handles concatenated gzip streams, such as Python `gzip`. Input order is retained; target positions are not globally sorted.

After an interrupted job has ended, rerun the same production command with the same profile. Committed chunk directories are validated by SHA256 and reused. Chunk directories become visible only after successful compression and an atomic directory rename. Old incomplete work directories are ignored. Completed traits are reused after final-file integrity checks.

**Ordinary gzip does not provide constant-time restart at a row offset:** the source prefix may be decompressed again to reach the saved uncompressed offset. Completed chunks are not remapped or regenerated. Changing input metadata, references, chain, executable or chunk size triggers a configuration mismatch; use a separate output directory for a different run.

Keep the chunk directories for resume and final accounting. This design uses disk space for both chunk and merged outputs. Progress percentages estimate compressed-input bytes and can reach 100% before the last compression/merge stage; completion is established by the manifests.

## SLURM

The original production requested one job with four CPUs and 24 GiB of memory. These are **requested resources**, not measured minimum requirements or a performance benchmark. The helper defaults to the original cluster's `bi2` partition; choose a partition valid on your cluster.

Submit:

```bash
GWAMA_PARTITION=bi2 bash hpc/submit.sh config.local.json
```

Check scheduler state **and data outputs** independently:

```bash
bash hpc/status.sh config.local.json
```

For direct commands with a returned job ID:

```bash
squeue -j JOB_ID -o "%.18i %.12P %.20j %.10T %.6C %.12M %.30R"
sacct -j JOB_ID --format=JobID,Partition,State,Elapsed,AllocCPUS,MaxRSS
python3 scripts/Step_03_liftover.py status --config config.local.json
```

## Recorded production outcomes

The following are the author's reported outcomes from the completed 14 November 2026 run; they are not benchmarks reproduced by the small bundled tests.

| Trait | Input records | Accepted SNPs | Review records | SNP retention | FASTA pair match |
|---|---:|---:|---:|---:|---:|
| MCH | 50,266,654 | 47,239,225 | 3,027,429 | 99.942900% | 99.998321% |
| MCV | 50,615,308 | 47,587,648 | 3,027,660 | 99.943068% | 99.998308% |
| RDW | 49,860,007 | 46,832,619 | 3,027,388 | 99.943369% | 99.998337% |

SNP retention is `PASS / ELIGIBLE_SNP`; FASTA pair match is `PASS / FASTA_CHECKED`. The 9,001,714 records labelled non-SNP/unsupported allele are still review records. Per-trait record totals overlap across traits and must not be interpreted as counts of distinct loci. All detailed counters, pilot outcomes and the reported successful final-audit status are retained in [results](results/README.md).

## Tests and reproducibility

Run offline standard-library tests:

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

Run the native UCSC integration regression on generated miniature genomes:

```bash
LIFTOVER_BIN="$(pwd)/tools/liftOver" python3 tests/native_integration.py
```

The native test covers plus/minus chains, ambiguous/unmapped records, FASTA mismatches, indels and malformed rows, unchanged statistics, gzip assembly, interruption/resume, corruption guards and locking. Human assembly-size preflight is independently tested; only that guard is bypassed for the miniature-genome integration fixture. This is not a reanalysis of the full three-trait data.

The **byte-identical scripts used for the reported production** are under [reproduction/original](reproduction/original). Their hashes are recorded in [reproduction/archive_sha256.json](reproduction/archive_sha256.json). The public portable scripts preserve the scientific rules but use explicit input mappings and configuration paths; they have a separate configuration version and should use a separate output directory. The public final auditor can also read the original completed run.

To export your own run metadata with absolute paths removed:

```bash
python3 scripts/export_provenance.py --root results/run --out results/provenance.redacted.json
```

Review the exported metadata before publishing. The original configuration ID still identifies the unredacted configuration; the redacted JSON is a provenance summary, not a replacement for the run's resume manifest.

## Manuscript and repository use

[docs/Supplementary_Methods_GWAMA_Liftover.docx](docs/Supplementary_Methods_GWAMA_Liftover.docx) contains the manuscript-ready methods text; an editable text version is in [docs/supplementary_methods.md](docs/supplementary_methods.md). Additional details are in the schema, references and reproduction notes.

To publish, extract the archive and upload the **contents of this repository folder** to a new GitHub repository, including `.gitignore` and `.github`. Alternatively, initialise a local Git repository and push to your chosen remote. Add the final repository URL to `CITATION.cff` after creating it; no repository URL or DOI is presumed here. Raw data and run outputs are excluded by `.gitignore`.

Original code and documentation are provided under the **MIT licence**. External data, UCSC resources and reference genomes retain their own terms. Cite the original data study, GWAMA and UCSC resources as applicable; see [docs/references.md](docs/references.md).
