#!/usr/bin/env bash
# Inspect local GWAS inputs and existing resources without changing data.
python3 -u - <<'PY'
from pathlib import Path
import bz2
import csv
import gzip
import lzma
import os
import platform
import shutil

base = Path.home() / "Drug_rankGWAS"
raw = base / "data/GWAS/raw"
print("=== Step 01 | GWAS liftover input audit ===")
print("Input:", raw)
print("Host:", platform.node(), "| System:", platform.system(), platform.machine())
print("CPU:", os.cpu_count(), "| Load:", os.getloadavg())
meminfo = Path("/proc/meminfo")
if meminfo.is_file():
    for line in meminfo.read_text().splitlines():
        if line.startswith(("MemTotal:", "MemAvailable:")):
            key, value, *_ = line.split()
            print(key, f"{int(value) / 1024**2:.1f} GiB")
if not raw.is_dir():
    raise SystemExit(f"[FAIL] Directory not found: {raw}")
print("Disk free:", f"{shutil.disk_usage(raw).free / 1024**3:.1f} GiB")

files = sorted(p for p in raw.iterdir() if p.is_file() and not p.name.startswith("."))
print(f"\n=== Directory files: {len(files)} ===")
for p in files:
    print(f"{p.stat().st_size / 1024**3:8.3f} GiB  {p.name}")

suffixes = (".gz", ".bgz", ".bgzip", ".bz2", ".xz", ".txt", ".tsv", ".csv", ".sumstats")
inputs = [p for p in files if p.name.lower().endswith(suffixes)]
print(f"\n=== Header previews: {len(inputs)} candidate files ===")
for p in inputs[:12]:
    print(f"\nFILE: {p.name}")
    try:
        with p.open("rb") as handle:
            magic = handle.read(6)
        opener = gzip.open if magic.startswith(b"\x1f\x8b") else (
            bz2.open if magic.startswith(b"BZh") else (
                lzma.open if magic.startswith(b"\xfd7zXZ\x00") else open))
        with opener(p, "rt", encoding="utf-8-sig", errors="replace") as handle:
            for _ in range(100):
                header = handle.readline().rstrip("\r\n")
                if header.strip() and not header.startswith("##"):
                    break
            else:
                raise ValueError("No header found within the first 100 lines")
            columns = (header.split("\t") if "\t" in header else
                       next(csv.reader([header])) if "," in header else header.split())
            print("N_COLUMNS:", len(columns))
            print("HEADER:", header)
            for i in range(1, 3):
                row = handle.readline().rstrip("\r\n")
                if row:
                    print(f"ROW_{i}:", row)
    except Exception as error:
        print(f"[READ_FAIL] {type(error).__name__}: {error}")
if len(inputs) > 12:
    print("[NOTE] Preview limited to the first 12 candidate files.")

print("\n=== Tools on PATH ===")
for name in ("python3", "liftOver", "samtools", "bgzip", "pigz"):
    print(f"{name}: {shutil.which(name) or 'NOT_ON_PATH'}")

print("\n=== Existing resource candidates ===")
atlas = Path.home() / "GRN_GWAS_Atlas"
resources = [
    base / "tools/liftOver",
    base / "reference/liftover/hg19ToHg38.over.chain.gz",
    atlas / "tools/ucsc_liftover_env/bin/liftOver",
    atlas / "reference/liftover/hg19ToHg38.over.chain.gz",
    atlas / "reference/GRCh38/GRCh38.primary_assembly.genome.fa",
    atlas / "reference/GRCh38/GRCh38.primary_assembly.genome.fa.fai",
]
for p in resources:
    print(f"[{'FOUND' if p.is_file() else 'MISSING'}] {p}")
print("\n[DONE] Preview only; no input files or running jobs were changed.")
print("[NOTE] A header preview does not verify the genome build or gzip integrity.")
PY
