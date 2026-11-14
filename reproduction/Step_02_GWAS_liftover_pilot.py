#!/usr/bin/env python3
"""Pilot coordinate and SNP allele checks; never modify the raw GWAS files."""
import csv
import gzip
import hashlib
import itertools
import json
import os
import re
import subprocess
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path.home() / "Drug_rankGWAS"
ATLAS = Path.home() / "GRN_GWAS_Atlas"
N_ROWS = 10000
PRIMARY = {"chr" + str(x) for x in range(1, 23)} | {"chrX", "chrY"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
STATUS = ("PASS", "NO_MAP", "MULTI_MAP", "NON_SNP_OR_INVALID_ALLELE", "BAD_ID",
          "BAD_SOURCE_COORD", "ID_ALLELE_MISMATCH", "NON_PRIMARY_OR_NON_1BP",
          "FASTA_N", "ALLELE_MISMATCH")


def validate_chain(path):
    with gzip.open(path, "rt", encoding="ascii") as handle:
        headers = [line.split() for line in handle if line.startswith("chain ")]
    if not any(len(h) == 13 and h[2:5] == ["chr1", "249250621", "+"]
               and h[7:9] == ["chr1", "248956422"] for h in headers):
        raise RuntimeError(f"Chain direction/assembly check failed: {path}")


def prepare_chain():
    path = BASE / "reference/liftover_ucsc/hg19ToHg38.over.chain.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        for host in ("hgdownload.soe.ucsc.edu", "hgdownload.cse.ucsc.edu"):
            url = f"https://{host}/goldenPath/hg19/liftOver/{path.name}"
            part = path.with_name(path.name + f".download.{os.getpid()}")
            print("[DOWNLOAD]", url, flush=True)
            try:
                with urllib.request.urlopen(url, timeout=30) as response, part.open("wb") as out:
                    total = 0
                    while True:
                        block = response.read(65536)
                        if not block:
                            break
                        total += len(block)
                        if total > 10000000:
                            raise RuntimeError("Unexpected chain download size")
                        out.write(block)
                validate_chain(part)
                os.replace(part, path)
                break
            except Exception as error:
                print(f"[DOWNLOAD_FAIL] {error}", flush=True)
            finally:
                part.unlink(missing_ok=True)
        if not path.exists():
            raise RuntimeError("Official chain download failed on both UCSC hosts")
    validate_chain(path)
    print("[CHAIN]", path, flush=True)
    print("[CHAIN_SHA256]", hashlib.sha256(path.read_bytes()).hexdigest(), flush=True)
    return path


class Fasta:
    def __init__(self, path):
        self.index = {}
        with Path(str(path) + ".fai").open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                self.index[fields[0]] = tuple(map(int, fields[1:5]))
        self.handle = path.open("rb")
        if self.handle.read(1) != b">":
            self.handle.close()
            raise RuntimeError("Expected an uncompressed FASTA beginning with '>'")
        for chrom, length in (("chr1", 248956422), ("chr10", 133797422)):
            key = chrom if chrom in self.index else chrom[3:]
            if key not in self.index or self.index[key][0] != length:
                self.handle.close()
                raise RuntimeError(f"GRCh38 FASTA chromosome-length check failed: {chrom}")

    def base(self, chrom, pos):
        key = chrom if chrom in self.index else chrom[3:]
        if key not in self.index:
            raise RuntimeError(f"FASTA contig missing: {chrom}")
        length, offset, bases, width = self.index[key]
        if not (1 <= pos <= length) or bases < 1 or width < bases:
            raise RuntimeError(f"Invalid FASTA coordinate/index: {chrom}:{pos}")
        position = pos - 1
        address = offset + (position // bases) * width + position % bases
        value = os.pread(self.handle.fileno(), 1, address).decode("ascii").upper()
        if value not in ("A", "C", "G", "T", "N"):
            raise RuntimeError(f"Invalid FASTA byte; check the .fai: {chrom}:{pos}")
        return value


def pilot(trait, lift, chain, fasta, output):
    source = BASE / f"data/GWAS/raw/{trait}_EA_GWAMA_20190214_GC.gz"
    bed = output / f"{trait}.hg19.bed"
    mapped = output / f"{trait}.hg38.bed"
    unmapped = output / f"{trait}.unmapped.bed"
    rows, eligible = [], set()
    print(f"\n[START] {trait}: first {N_ROWS:,} data rows", flush=True)
    with gzip.open(source, "rt", encoding="utf-8-sig") as handle, bed.open("w") as out:
        header = handle.readline().rstrip("\r\n")
        sep = "\t" if "\t" in header else None
        names = header.split(sep)
        required = {"rs_number", "reference_allele", "other_allele", "beta", "eaf"}
        if len(names) != 19 or not required.issubset(names) or len(set(names)) != 19:
            raise RuntimeError(f"Unexpected GWAS header: {source}")
        for rid, line in enumerate(itertools.islice(handle, N_ROWS), 1):
            fields = line.rstrip("\r\n").split(sep)
            if len(fields) != len(names):
                raise RuntimeError(f"Column-count mismatch: {trait}, data row {rid}")
            record = dict(zip(names, fields))
            extra = {"source_row": rid, "qc": "BAD_ID"}
            match = re.fullmatch(r"(?:chr)?([0-9]+|X|Y):([0-9]+)_([^_]+)_([^_]+)", record["rs_number"])
            if match:
                chrom, pos, a, b = match.groups()
                chrom, pos = "chr" + chrom, int(pos)
                a, b = a.upper(), b.upper()
                ra, oa = record["reference_allele"].upper(), record["other_allele"].upper()
                extra.update(chr37=chrom, pos37=pos)
                if chrom not in PRIMARY or pos < 1:
                    extra["qc"] = "BAD_SOURCE_COORD"
                elif any(x not in ("A", "C", "G", "T") for x in (a, b, ra, oa)) or a == b:
                    extra["qc"] = "NON_SNP_OR_INVALID_ALLELE"
                elif {a, b} != {ra, oa}:
                    extra["qc"] = "ID_ALLELE_MISMATCH"
                else:
                    extra["qc"] = "NO_MAP"
                    eligible.add(rid)
                    out.write(f"{chrom}\t{pos-1}\t{pos}\t{rid}\t0\t+\n")
            rows.append((record, extra))
    if not rows:
        raise RuntimeError(f"No data rows: {source}")
    mappings = defaultdict(set)
    if eligible:
        command = [str(lift), "-minMatch=1.0", "-multiple", str(bed), str(chain), str(mapped), str(unmapped)]
        with (output / f"{trait}.liftOver.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)
        with mapped.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) != 6 or fields[5] not in ("+", "-"):
                    raise RuntimeError(f"Unexpected mapped BED6 row: {line!r}")
                chrom, start, end, rid, _, strand = fields
                rid = int(rid)
                if rid not in eligible:
                    raise RuntimeError(f"Unknown mapped row ID: {rid}")
                mappings[rid].add((chrom, int(start), int(end), strand))
    counts = Counter()
    for record, extra in rows:
        destinations = mappings.get(extra["source_row"], set())
        if len(destinations) > 1:
            extra["qc"] = "MULTI_MAP"
        elif len(destinations) == 1:
            chrom, start, end, strand = next(iter(destinations))
            extra.update(chr38=chrom, pos38=start + 1, lift_strand=strand)
            if chrom not in PRIMARY or end - start != 1:
                extra["qc"] = "NON_PRIMARY_OR_NON_1BP"
            else:
                ra, oa = record["reference_allele"].upper(), record["other_allele"].upper()
                if strand == "-":
                    ra, oa = ra.translate(COMPLEMENT), oa.translate(COMPLEMENT)
                base = fasta.base(chrom, start + 1)
                role = "reference_allele" if base == ra else "other_allele" if base == oa else "NEITHER"
                extra.update(reference_allele_hg38=ra, other_allele_hg38=oa,
                             fasta_base_hg38=base, fasta_matches_column=role)
                extra["qc"] = "PASS" if role != "NEITHER" else "FASTA_N" if base == "N" else "ALLELE_MISMATCH"
                counts["FASTA_CHECKED"] += 1
                counts["MINUS_CHECKED"] += strand == "-"
                counts["PLUS_CHECKED"] += strand == "+"
                counts["MINUS_PASS"] += strand == "-" and extra["qc"] == "PASS"
                counts["PLUS_PASS"] += strand == "+" and extra["qc"] == "PASS"
        counts[extra["qc"]] += 1
    columns = ["source_row", "chr37", "pos37", "chr38", "pos38", "lift_strand",
               "reference_allele_hg38", "other_allele_hg38", "fasta_base_hg38", "fasta_matches_column", "qc"]
    with (output / f"{trait}.pilot.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names + columns, delimiter="\t", restval="NA")
        writer.writeheader()
        writer.writerows({**record, **extra} for record, extra in rows)
    result = {"trait": trait, "SAMPLED": len(rows), "ELIGIBLE_SNP": len(eligible)}
    result.update({key: counts[key] for key in STATUS + ("FASTA_CHECKED", "PLUS_CHECKED", "MINUS_CHECKED")})
    result["PAIR_MATCH_RATE"] = f"{counts['PASS'] / counts['FASTA_CHECKED']:.6f}" if counts["FASTA_CHECKED"] else "NA"
    for strand in ("PLUS", "MINUS"):
        result[strand + "_MATCH_RATE"] = (f"{counts[strand + '_PASS'] / counts[strand + '_CHECKED']:.6f}"
                                         if counts[strand + "_CHECKED"] else "NA")
    print(json.dumps(result, indent=2), flush=True)
    return result


def main():
    lift = ATLAS / "tools/ucsc_liftover_env/bin/liftOver"
    fasta_path = ATLAS / "reference/GRCh38/GRCh38.primary_assembly.genome.fa"
    if not lift.is_file() or not os.access(lift, os.X_OK):
        raise RuntimeError(f"liftOver is missing or not executable: {lift}")
    usage = subprocess.run([str(lift)], capture_output=True, text=True, timeout=15)
    if "usage:" not in (usage.stdout + usage.stderr).lower():
        raise RuntimeError(f"liftOver launch failed: {usage.stdout}{usage.stderr}")
    parent = BASE / "processed/GWAS/liftover_pilot"
    parent.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="Step_02_", dir=parent))
    print("[OUTPUT]", output, flush=True)
    (output / "liftOver.usage.txt").write_text(usage.stdout + usage.stderr)
    chain = prepare_chain()
    fasta = Fasta(fasta_path)
    try:
        results = [pilot(trait, lift, chain, fasta, output) for trait in ("MCH", "MCV", "RDW")]
    finally:
        fasta.handle.close()
    names = list(dict.fromkeys(key for item in results for key in item))
    with (output / "QC_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", restval=0)
        writer.writeheader()
        writer.writerows(results)
    (output / "resources.json").write_text(json.dumps({"sample_rows": N_ROWS, "liftOver": str(lift),
        "chain": str(chain), "chain_sha256": hashlib.sha256(chain.read_bytes()).hexdigest(),
        "GRCh38_fasta": str(fasta_path), "scope": "First rows only; SNP allele checks; not full conversion"}, indent=2))
    print("\n[PILOT_COMPLETE]", output, flush=True)
    print("[QC]", output / "QC_summary.tsv", flush=True)
    print("[NOTE] Review QC before full conversion. Raw data and statistical columns were not changed.", flush=True)


if __name__ == "__main__":
    main()
