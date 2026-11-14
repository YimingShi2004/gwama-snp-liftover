#!/usr/bin/env python3
"""Verify completed SNP liftover files, provenance fingerprints and row accounting.

This audit checks file integrity and production accounting. Genome-wide target
duplicates, GWAS statistical QC and LD-panel harmonization are separate steps.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
import tempfile
import time
import traceback
from collections import Counter
from pathlib import Path

TRAITS = ("MCH", "MCV", "RDW")
STATUSES = ("PASS", "BAD_COLUMN_COUNT", "BAD_ID", "BAD_SOURCE_COORD",
            "NON_SNP_NEEDS_VARIANT_LIFTOVER", "ID_ALLELE_MISMATCH", "NO_MAP",
            "MULTI_MAP", "NON_PRIMARY_OR_NON_1BP", "FASTA_N", "ALLELE_MISMATCH")
PRE_MAP = STATUSES[1:6]
POST_MAP = STATUSES[6:]


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def load(path):
    with path.open() as handle:
        return json.load(handle)


def sha256(path, verbose=False):
    digest = hashlib.sha256()
    done, total, last = 0, path.stat().st_size, time.monotonic()
    if verbose:
        print(f"[HASH] {path.name}: {total / 1024**3:.3f} GiB", flush=True)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b""):
            digest.update(block)
            done += len(block)
            if verbose and time.monotonic() - last >= 5:
                print(f"  {path.name}: {100 * done / max(total, 1):.1f}%", flush=True)
                last = time.monotonic()
    return digest.hexdigest()


def check_signature(saved):
    path = Path(saved["path"])
    check(path.is_file(), f"Provenance file missing: {path}")
    st = path.stat()
    check(st.st_size == saved["size"] and st.st_mtime_ns == saved["mtime_ns"],
          f"Provenance file size/mtime changed: {path}")
    if "sha256" in saved:
        check(sha256(path) == saved["sha256"], f"Provenance SHA256 mismatch: {path}")


def same_counts(first, second, label):
    check(all(first.get(k, 0) == second.get(k, 0) for k in set(first) | set(second)),
          f"Counter mismatch: {label}")


def check_counts(counts, label):
    c = Counter(counts)
    check(all(isinstance(v, int) and v >= 0 for v in counts.values()), f"Invalid count: {label}")
    checks = {
        "INPUT=PASS+REVIEW": c["INPUT"] == c["PASS"] + c["REVIEW"],
        "INPUT=sum(statuses)": c["INPUT"] == sum(c[k] for k in STATUSES),
        "INPUT=ELIGIBLE+pre-map exclusions": c["INPUT"] == c["ELIGIBLE_SNP"] + sum(c[k] for k in PRE_MAP),
        "ELIGIBLE=PASS+post-map exclusions": c["ELIGIBLE_SNP"] == c["PASS"] + sum(c[k] for k in POST_MAP),
        "FASTA=PASS+N+mismatch": c["FASTA_CHECKED"] == c["PASS"] + c["FASTA_N"] + c["ALLELE_MISMATCH"],
        "FASTA=PLUS+MINUS": c["FASTA_CHECKED"] == c["PLUS_CHECKED"] + c["MINUS_CHECKED"],
        "PASS=PLUS_PASS+MINUS_PASS": c["PASS"] == c["PLUS_PASS"] + c["MINUS_PASS"],
        "strand pass bounds": c["PLUS_PASS"] <= c["PLUS_CHECKED"] and c["MINUS_PASS"] <= c["MINUS_CHECKED"],
        "palindromic pass bounds": c["PALINDROMIC_PASS"] <= c["PASS"],
    }
    for name, passed in checks.items():
        check(passed, f"Row-accounting failure in {label}: {name}")
    return c


def check_summary_row(row, counts, trait):
    for key, value in row.items():
        if key in ("PAIR_MATCH_RATE", "PLUS_MATCH_RATE", "MINUS_MATCH_RATE"):
            prefix = key.split("_")[0]
            n = counts["FASTA_CHECKED"] if prefix == "PAIR" else counts[prefix + "_CHECKED"]
            p = counts["PASS"] if prefix == "PAIR" else counts[prefix + "_PASS"]
            check(value == "NA" if not n else abs(float(value) - p / n) <= 5.1e-9,
                  f"QC summary rate mismatch: {trait}, {key}")
        elif key != "trait":
            check(int(value) == counts[key], f"QC summary count mismatch: {trait}, {key}")


def inspect_header(path, accepted):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"rs_number_GRCh37", "reference_allele_GRCh37", "other_allele_GRCh37", "beta", "se", "eaf",
                    "p-value", "chr38", "pos38", "lift_strand", "reference_allele_hg38", "other_allele_hg38",
                    "REF_hg38", "ALT_hg38", "variant_id_hg38", "ALT_hg38_from", "liftover_qc", "source_row"}
        check(required.issubset(fields) and len(fields) == len(set(fields)), f"Unexpected header: {path}")
        if not accepted:
            check("raw_line_base64" in fields, f"Review payload column missing: {path}")
        first = next(reader, None)
        if accepted and first:
            check(None not in first and first["liftover_qc"] == "PASS", f"Malformed first accepted row: {path}")
            check(first["lift_strand"] in ("+", "-") and int(first["pos38"]) >= 1, f"Invalid coordinate: {path}")
            ref, alt = first["REF_hg38"], first["ALT_hg38"]
            check(ref in ("A", "C", "G", "T") and alt in ("A", "C", "G", "T") and ref != alt,
                  f"Invalid first SNP: {path}")
            check({ref, alt} == {first["reference_allele_hg38"], first["other_allele_hg38"]},
                  f"First allele pair inconsistent: {path}")
            check(first["variant_id_hg38"] == f"{first['chr38']}:{first['pos38']}:{ref}:{alt}",
                  f"First variant ID inconsistent: {path}")
    return len(fields), first is not None


def audit(root):
    check(root.is_dir(), f"Output directory missing: {root}")
    print("[ROOT]", root, flush=True)
    manifest, overall = load(root / "run_manifest.json"), load(root / "COMPLETE.json")
    configuration = manifest["configuration"]
    traits = tuple(configuration.get("traits", overall["traits"]))
    check(bool(traits), "No traits recorded in manifest")
    config_id = hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).hexdigest()
    check(config_id == manifest["config_id"] == overall["config_id"], "Run configuration ID mismatch")
    check(configuration["version"] in ("Drug_rankGWAS_UCSC_SNP_v1", "GWAMA_SNP_liftover_v1"), "Unexpected production version")
    print("[CHECK] Input/resource fingerprints", flush=True)
    for saved in configuration["source_files"]:
        check_signature(saved)
    for key in ("chain", "liftOver", "fasta", "fai"):
        check_signature(configuration[key])
    with (root / "QC_summary.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    check(len(rows) == len(traits) and {r["trait"] for r in rows} == set(traits), "QC summary trait list mismatch")
    summary_rows = {r["trait"]: r for r in rows}
    results = []
    for trait in traits:
        directory = root / trait
        complete = load(directory / "COMPLETE.json")
        check(complete["config_id"] == config_id and complete["trait"] == trait, f"Manifest mismatch: {trait}")
        check(complete["input_gzip_crc_checked_to_eof"] is True, f"Production did not record source EOF/CRC check: {trait}")
        counts = check_counts(complete["counts"], trait)
        same_counts(counts, overall["traits"][trait], trait + " root/trait manifests")
        chunks = sorted(p for p in (directory / "chunks").iterdir() if re.fullmatch(r"chunk_[0-9]{6}", p.name))
        check(len(chunks) == complete["chunks"] and bool(chunks), f"Checkpoint count mismatch: {trait}")
        accumulated, end, gzip_offset = Counter(), 0, 0
        for i, path in enumerate(chunks):
            check(path.name == f"chunk_{i:06d}", f"Checkpoint sequence gap: {path}")
            q = load(path / "qc.json")
            check(q["config_id"] == config_id and q["row_start"] == end + 1, f"Checkpoint provenance/sequence mismatch: {path}")
            current = check_counts(q["counts"], path.name)
            check(current["INPUT"] == q["row_end"] - q["row_start"] + 1, f"Checkpoint row range mismatch: {path}")
            check(q["gzip_uncompressed_end"] > gzip_offset, f"Invalid source offset sequence: {path}")
            accumulated.update(current)
            end, gzip_offset = q["row_end"], q["gzip_uncompressed_end"]
        same_counts(accumulated, counts, trait + " checkpoint totals")
        check(end == counts["INPUT"], f"Final input row mismatch: {trait}")
        check_summary_row(summary_rows[trait], counts, trait)
        expected = {f"{trait}.SNP.GRCh38.tsv.gz", f"{trait}.review.tsv.gz"}
        check(set(complete["outputs"]) == expected, f"Unexpected final file list: {trait}")
        for name in sorted(expected):
            path, saved = directory / name, complete["outputs"][name]
            check(path.is_file() and path.stat().st_size == saved["bytes"], f"Final file size mismatch: {path}")
            check(sha256(path, verbose=True) == saved["sha256"], f"Final file SHA256 mismatch: {path}")
            accepted = name.endswith(".SNP.GRCh38.tsv.gz")
            nfields, has_data = inspect_header(path, accepted)
            expected_rows = counts["PASS"] if accepted else counts["REVIEW"]
            check(has_data == (expected_rows > 0), f"Final file empty/nonempty mismatch: {path}")
            print(f"[FILE_PASS] {name}: SHA256 matched, columns={nfields}", flush=True)
        print(f"[TRAIT_PASS] {trait}: INPUT={counts['INPUT']:,} = PASS {counts['PASS']:,} + REVIEW {counts['REVIEW']:,}", flush=True)
        result = dict(summary_rows[trait])
        result.update(audit="PASS", final_files_sha256="PASS", row_accounting="PASS")
        results.append(result)
    report = {"audit": "SNP_LIFTOVER_FINAL_AUDIT_PASS", "audited_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "config_id": config_id, "results": results,
              "scope": "Final compressed-file SHA256/size, fingerprints, checkpoint accounting and output schema; source CRC verified during production",
              "next_stage": "Genome-wide duplicate checks, GWAS statistical QC, LD harmonization and GRN input compatibility",
              "non_snp_scope": "Non-SNP/unsupported-allele rows remain in review and are not certified GRCh38 variants"}
    parent = root / "final_audit"
    parent.mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="Step_04_", dir=parent))
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    with (output / "audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    print("\n[SNP_LIFTOVER_FINAL_AUDIT_PASS]", flush=True)
    print("[REPORT]", output / "audit.tsv", flush=True)
    print("[SCOPE] Validated SNP liftover files; review rows retained for separate resolution.", flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()
    if args.root is None:
        if args.config:
            profile = args.config.expanduser().resolve()
            value = load(profile).get("out", "results/run")
            args.root = Path(value).expanduser()
            if not args.root.is_absolute():
                args.root = profile.parent / args.root
        else:
            args.root = Path.cwd() / "results/run"
    try:
        audit(args.root.expanduser().resolve())
    except Exception:
        print("[FINAL_AUDIT_FAILED] See the error below; existing results were preserved.", file=sys.stderr, flush=True)
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
