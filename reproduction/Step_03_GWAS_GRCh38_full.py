#!/usr/bin/env python3
"""Stream three GWAS files through UCSC liftOver with SNP QC and atomic checkpoints.

Accepted files contain validated biallelic SNPs. Every other input row is kept
in a review file. Non-SNPs require a separate variant-aware conversion.
Only allele complementation is performed; statistical values are never changed.
"""
import argparse
import base64
import csv
import fcntl
import gzip
import hashlib
import itertools
import json
import mmap
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from collections import Counter
from pathlib import Path

VERSION = "Drug_rankGWAS_UCSC_SNP_v1"
TRAITS = ("MCH", "MCV", "RDW")
PRIMARY = {"chr" + str(c) for c in range(1, 23)} | {"chrX", "chrY"}
PATTERN = re.compile(rb"(?:chr)?([0-9]+|X|Y):([0-9]+)_([^_]+)_([^_]+)")
TRANSLATE = bytes.maketrans(b"ACGT", b"TGCA")
STATUS = ("PASS", "BAD_COLUMN_COUNT", "BAD_ID", "BAD_SOURCE_COORD",
          "NON_SNP_NEEDS_VARIANT_LIFTOVER", "ID_ALLELE_MISMATCH", "NO_MAP",
          "MULTI_MAP", "NON_PRIMARY_OR_NON_1BP", "FASTA_N", "ALLELE_MISMATCH")
EXTRA = ("source_row", "chr37", "pos37", "chr38", "pos38", "lift_strand",
         "reference_allele_hg38", "other_allele_hg38", "REF_hg38", "ALT_hg38",
         "variant_id_hg38", "ALT_hg38_from", "liftover_qc", "mapping_candidates")
RENAME = {"rs_number": "rs_number_GRCh37", "reference_allele": "reference_allele_GRCh37",
          "other_allele": "other_allele_GRCh37"}


def stamp():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def log(message):
    print(f"[{stamp()}] {message}", flush=True)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def signature(path):
    st = path.stat()
    return {"path": str(path.resolve()), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def atomic_json(path, value):
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def load_json(path):
    with path.open() as handle:
        return json.load(handle)


def file_info(path):
    return {"bytes": path.stat().st_size, "sha256": sha256(path)}


def verify(path, info):
    if not path.is_file() or path.stat().st_size != info["bytes"] or sha256(path) != info["sha256"]:
        raise RuntimeError(f"Checkpoint integrity check failed: {path}; preserve it for diagnosis")


class Progress:
    def __init__(self, root):
        self.path = root / "progress.json"
        self.state = {"version": VERSION, "pid": os.getpid(), "host": os.uname().nodename,
                      "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""), "started_utc": stamp()}

    def update(self, stage, **values):
        self.state.update(values, stage=stage, updated_utc=stamp())
        atomic_json(self.path, self.state)
        log(f"stage={stage} trait={self.state.get('trait', '-')} "
            f"committed_rows={self.state.get('committed_rows', 0):,}")


class Fasta:
    def __init__(self, path):
        self.index = {}
        with Path(str(path) + ".fai").open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                if len(fields) != 5 or fields[0] in self.index:
                    raise RuntimeError("Malformed or duplicated FASTA index entry")
                self.index[fields[0]] = tuple(map(int, fields[1:]))
        self.handle = path.open("rb")
        if self.handle.read(1) != b">":
            self.handle.close()
            raise RuntimeError("Expected an uncompressed FASTA")
        self.mapping = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)

    def key(self, chrom):
        return chrom if chrom in self.index else chrom[3:]

    def check_length(self, chrom, length):
        key = self.key(chrom)
        if key not in self.index or self.index[key][0] != length:
            raise RuntimeError(f"FASTA/chain chromosome length mismatch: {chrom}")

    def base(self, chrom, position):
        key = self.key(chrom)
        if key not in self.index:
            raise RuntimeError(f"Missing FASTA chromosome: {chrom}")
        length, offset, bases, width = self.index[key]
        if not 1 <= position <= length or bases < 1 or width < bases:
            raise RuntimeError(f"Invalid FASTA coordinate/index: {chrom}:{position}")
        pos = position - 1
        address = offset + (pos // bases) * width + pos % bases
        if not 0 <= address < len(self.mapping):
            raise RuntimeError(f"FASTA index points outside the file: {chrom}:{position}")
        base = bytes((self.mapping[address],)).upper()
        if base not in (b"A", b"C", b"G", b"T", b"N"):
            raise RuntimeError(f"Invalid FASTA byte at {chrom}:{position}; check the .fai")
        return base

    def close(self):
        self.mapping.close()
        self.handle.close()


def inspect_chain(path, fasta):
    source_sizes = {}
    target_sizes = {}
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            if not line.startswith("chain "):
                continue
            h = line.split()
            if len(h) != 13 or h[4] != "+":
                raise RuntimeError("Unsupported or malformed chain header")
            for sizes, name, length in ((source_sizes, h[2], int(h[3])), (target_sizes, h[7], int(h[8]))):
                if name in sizes and sizes[name] != length:
                    raise RuntimeError("Inconsistent chain chromosome length")
                sizes[name] = length
    if source_sizes.get("chr1") != 249250621 or target_sizes.get("chr1") != 248956422:
        raise RuntimeError("Expected a direct hg19 -> hg38 chain")
    if source_sizes.get("chr10") != 135534747 or target_sizes.get("chr10") != 133797422:
        raise RuntimeError("Chain direction/assembly check failed on chr10")
    for chrom in PRIMARY:
        if chrom not in source_sizes or chrom not in target_sizes:
            raise RuntimeError(f"Primary chromosome missing from chain: {chrom}")
        fasta.check_length(chrom, target_sizes[chrom])
    return source_sizes


def run_logged(command, logfile, progress, stage):
    progress.update(stage)
    with logfile.open("ab") as handle:
        handle.write(("COMMAND: " + shlex.join(map(str, command)) + "\n").encode())
        handle.flush()
        process = subprocess.Popen(list(map(str, command)), stdout=handle, stderr=handle)
        try:
            while True:
                try:
                    code = process.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    progress.update(stage + "_running")
            if code:
                raise RuntimeError(f"Command failed ({code}); see {logfile}")
        except BaseException:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise


def compress_file(source, destination, pigz, threads, progress):
    if pigz:
        progress.update("compress_" + source.stem)
        with source.open("rb") as inp, destination.open("wb") as out:
            proc = subprocess.Popen([pigz, "-p", str(threads), "-1", "-c"], stdin=inp, stdout=out)
            try:
                while True:
                    try:
                        code = proc.wait(timeout=30)
                        break
                    except subprocess.TimeoutExpired:
                        progress.update("compress_running")
                if code:
                    raise RuntimeError(f"pigz failed ({code}): {source}")
            except BaseException:
                if proc.poll() is None:
                    proc.terminate()
                    proc.wait(timeout=10)
                raise
    else:
        progress.update("compress_python_gzip")
        with source.open("rb") as inp, destination.open("wb") as out:
            with gzip.GzipFile(filename="", mode="wb", compresslevel=1, fileobj=out, mtime=0) as gz:
                shutil.copyfileobj(inp, gz, length=4 * 1024**2)
    source.unlink()
    return file_info(destination)


def parse_record(raw, rowid, columns, sep, source_sizes):
    fields = raw.split(sep)
    if len(fields) != len(columns):
        return (raw, rowid, "NA", 0, b"", b"", "BAD_COLUMN_COUNT")
    marker = fields[columns["rs_number"]]
    match = PATTERN.fullmatch(marker)
    if not match:
        return (raw, rowid, "NA", 0, b"", b"", "BAD_ID")
    chrom, position, a, b = match.groups()
    chrom = "chr" + (str(int(chrom)) if chrom.isdigit() else chrom.decode("ascii"))
    position = int(position)
    ra = fields[columns["reference_allele"]].upper()
    oa = fields[columns["other_allele"]].upper()
    if chrom not in PRIMARY or not 1 <= position <= source_sizes.get(chrom, 0):
        status = "BAD_SOURCE_COORD"
    elif any(x not in (b"A", b"C", b"G", b"T") for x in (a.upper(), b.upper(), ra, oa)) or a.upper() == b.upper():
        status = "NON_SNP_NEEDS_VARIANT_LIFTOVER"
    elif {a.upper(), b.upper()} != {ra, oa}:
        status = "ID_ALLELE_MISMATCH"
    else:
        status = "NO_MAP"
    return (raw, rowid, chrom, position, ra, oa, status)


def map_chunk(records, work, args, chain, progress):
    bed, mapped, missing = (work / name for name in ("input.bed", "mapped.bed", "unmapped.bed"))
    eligible = set()
    with bed.open("w") as handle:
        for _, rowid, chrom, pos, _, _, status in records:
            if status == "NO_MAP":
                eligible.add(rowid)
                handle.write(f"{chrom}\t{pos-1}\t{pos}\t{rowid}\t0\t+\n")
    destinations = {}
    if eligible:
        command = [args.lift, "-minMatch=1.0", "-multiple", bed, chain, mapped, missing]
        run_logged(command, work / "liftOver.log", progress, "liftOver")
        with mapped.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) != 6 or fields[5] not in ("+", "-"):
                    raise RuntimeError(f"Unexpected mapped BED6 output: {line!r}")
                chrom, start, end, rowid, _, strand = fields
                rowid = int(rowid)
                if rowid not in eligible:
                    raise RuntimeError(f"Unknown mapped row ID: {rowid}")
                value = (chrom, int(start), int(end), strand)
                previous = destinations.get(rowid)
                if previous is None:
                    destinations[rowid] = value
                elif isinstance(previous, tuple):
                    if previous != value:
                        destinations[rowid] = {previous, value}
                else:
                    previous.add(value)
    return destinations, len(eligible)


def write_rows(records, destinations, fasta, work, sep, ncolumns, progress):
    counts = Counter()
    accepted, review = work / "accepted.tsv", work / "review.tsv"
    last = time.monotonic()
    progress.update("FASTA_and_write")
    with accepted.open("wb", buffering=1024**2) as good, review.open("wb", buffering=1024**2) as bad:
        for i, record in enumerate(records, 1):
            raw, rowid, chrom37, pos37, ra, oa, status = record
            values = {key: "NA" for key in EXTRA}
            values.update(source_row=rowid, chr37=chrom37, pos37=pos37 or "NA")
            destination = destinations.get(rowid)
            if isinstance(destination, set):
                status = "MULTI_MAP"
                values["mapping_candidates"] = ";".join(f"{c}:{s+1}:{e}:{st}" for c, s, e, st in sorted(destination))
            elif destination is not None:
                chrom, start, end, strand = destination
                values.update(chr38=chrom, pos38=start + 1, lift_strand=strand)
                values["mapping_candidates"] = f"{chrom}:{start+1}:{end}:{strand}"
                if chrom not in PRIMARY or end - start != 1:
                    status = "NON_PRIMARY_OR_NON_1BP"
                else:
                    if strand == "-":
                        ra, oa = ra.translate(TRANSLATE), oa.translate(TRANSLATE)
                    base = fasta.base(chrom, start + 1)
                    values.update(reference_allele_hg38=ra.decode(), other_allele_hg38=oa.decode(), REF_hg38=base.decode())
                    counts["FASTA_CHECKED"] += 1
                    group = "PLUS" if strand == "+" else "MINUS"
                    counts[group + "_CHECKED"] += 1
                    if base in (ra, oa):
                        alt = oa if base == ra else ra
                        origin = RENAME["other_allele"] if base == ra else RENAME["reference_allele"]
                        values.update(ALT_hg38=alt.decode(), ALT_hg38_from=origin,
                                      variant_id_hg38=f"{chrom}:{start+1}:{base.decode()}:{alt.decode()}")
                        status = "PASS"
                        counts[group + "_PASS"] += 1
                        counts["PALINDROMIC_PASS"] += {ra, oa} in ({b"A", b"T"}, {b"C", b"G"})
                    else:
                        status = "FASTA_N" if base == b"N" else "ALLELE_MISMATCH"
            values["liftover_qc"] = status
            counts[status] += 1
            extra = b"\t".join(str(values[key]).encode("ascii") for key in EXTRA)
            malformed = status == "BAD_COLUMN_COUNT"
            body = b"\t".join([b"NA"] * ncolumns) if malformed else raw if sep == b"\t" else b"\t".join(raw.split())
            line = body + b"\t" + extra
            if status == "PASS":
                good.write(line + b"\n")
            else:
                payload = base64.b64encode(raw) if malformed else b"NA"
                bad.write(line + b"\t" + payload + b"\n")
            if i % 50000 == 0 and time.monotonic() - last >= 20:
                progress.update("FASTA_and_write", rows_in_current_chunk=i)
                last = time.monotonic()
    counts["INPUT"] = len(records)
    if sum(counts[key] for key in STATUS) != counts["INPUT"]:
        raise RuntimeError("Input/output row conservation check failed")
    counts["REVIEW"] = counts["INPUT"] - counts["PASS"]
    return counts, accepted, review


def header_info(handle):
    raw = handle.readline().rstrip(b"\r\n")
    raw = raw.removeprefix(b"\xef\xbb\xbf")
    sep = b"\t" if b"\t" in raw else None
    names = raw.decode("utf-8").split("\t" if sep else None)
    required = {"rs_number", "reference_allele", "other_allele", "beta", "se", "eaf", "p-value"}
    if len(names) != 19 or len(set(names)) != 19 or not required.issubset(names):
        raise RuntimeError("Input header differs from the audited 19-column GWAS format")
    output = [RENAME.get(name, name) for name in names] + list(EXTRA)
    return names, {name: i for i, name in enumerate(names)}, sep, output


def merged_file(path, header, chunks, key, progress):
    temporary = path.with_name(path.name + ".partial")
    digest = hashlib.sha256()
    log(f"MERGE {path.name}")
    progress.update("merge_" + key)
    last = time.monotonic()
    with temporary.open("wb") as handle:
        member = gzip.compress(("\t".join(header) + "\n").encode(), compresslevel=1, mtime=0)
        handle.write(member)
        digest.update(member)
        for chunk in chunks:
            with (chunk / (key + ".tsv.gz")).open("rb") as inp:
                for block in iter(lambda: inp.read(4 * 1024**2), b""):
                    handle.write(block)
                    digest.update(block)
                    if time.monotonic() - last >= 30:
                        progress.update("merge_" + key + "_running")
                        last = time.monotonic()
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def process_trait(trait, args, root, config_id, source_sizes, fasta, progress):
    source = args.base / f"data/GWAS/raw/{trait}_EA_GWAMA_20190214_GC.gz"
    directory = root / trait
    directory.mkdir(exist_ok=True)
    chunks_root = directory / "chunks"
    chunks_root.mkdir(exist_ok=True)
    done = directory / "COMPLETE.json"
    progress.update("trait_start", trait=trait, committed_rows=0, gzip_percent=0, counts={}, rows_in_current_chunk=0)
    if done.exists():
        saved = load_json(done)
        if saved["config_id"] != config_id:
            raise RuntimeError(f"Completion manifest configuration mismatch: {trait}")
        progress.update("verify_completed_outputs")
        for name, info in saved["outputs"].items():
            verify(directory / name, info)
        log(f"REUSE_FINAL {trait}: {saved['counts']['INPUT']:,} input rows")
        return saved
    completed = sorted(p for p in chunks_root.iterdir() if re.fullmatch(r"chunk_[0-9]{6}", p.name))
    counts = Counter()
    offset = 0
    for index, chunk in enumerate(completed):
        if chunk.name != f"chunk_{index:06d}":
            raise RuntimeError(f"Non-contiguous checkpoint sequence: {chunk}")
        qc = load_json(chunk / "qc.json")
        if (qc["config_id"] != config_id or qc["row_start"] != offset + 1
                or qc["row_end"] < qc["row_start"]
                or qc["counts"]["INPUT"] != qc["row_end"] - qc["row_start"] + 1
                or qc["counts"].get("PASS", 0) + qc["counts"].get("REVIEW", 0) != qc["counts"]["INPUT"]):
            raise RuntimeError(f"Checkpoint configuration/row sequence mismatch: {chunk}")
        progress.update("verify_checkpoint", chunk=index, committed_rows=offset)
        for name, info in qc["outputs"].items():
            verify(chunk / name, info)
        counts.update(qc["counts"])
        offset = qc["row_end"]
    source_signature = signature(source)
    with source.open("rb") as compressed, gzip.GzipFile(fileobj=compressed, mode="rb") as handle:
        names, columns, sep, output_header = header_info(handle)
        if completed:
            checkpoint = load_json(completed[-1] / "qc.json")
            log(f"RESUME {trait}: reuse {len(completed)} chunks / {offset:,} rows; decompress gzip prefix to checkpoint")
            progress.update("resume_decompress_prefix", committed_rows=offset)
            target = checkpoint["gzip_uncompressed_end"]
            if target < handle.tell():
                raise RuntimeError("Checkpoint offset precedes the GWAS header")
            last = time.monotonic()
            while handle.tell() < target:
                block = handle.read(min(8 * 1024**2, target - handle.tell()))
                if not block:
                    raise RuntimeError("Input ended before checkpoint offset")
                if time.monotonic() - last >= 30:
                    progress.update("resume_decompress_prefix", gzip_percent=100 * compressed.tell() / source_signature["size"])
                    last = time.monotonic()
        chunk_index = len(completed)
        while True:
            started = time.monotonic()
            progress.update("read_chunk", chunk=chunk_index, committed_rows=offset, rows_in_current_chunk=0)
            records = []
            last = time.monotonic()
            for i, raw in enumerate(itertools.islice(handle, args.chunk_size), 1):
                raw = raw.rstrip(b"\r\n")
                records.append(parse_record(raw, offset + i, columns, sep, source_sizes))
                if i % 50000 == 0 and time.monotonic() - last >= 20:
                    progress.update("read_chunk", rows_in_current_chunk=i)
                    last = time.monotonic()
            if not records:
                break
            end_offset = handle.tell()
            if shutil.disk_usage(root).free < args.min_free_gib * 1024**3:
                raise RuntimeError(f"Less than {args.min_free_gib:g} GiB free; existing checkpoints are preserved")
            work = Path(tempfile.mkdtemp(prefix=f".work_{chunk_index:06d}_", dir=directory))
            destinations, eligible = map_chunk(records, work, args, args.chain, progress)
            current, accepted, review = write_rows(records, destinations, fasta, work, sep, len(names), progress)
            current["ELIGIBLE_SNP"] = eligible
            outputs = {}
            for plain in (accepted, review):
                gz = plain.with_suffix(plain.suffix + ".gz")
                outputs[gz.name] = compress_file(plain, gz, args.pigz, args.compress_threads, progress)
            qc = {"version": VERSION, "config_id": config_id, "row_start": offset + 1,
                  "row_end": offset + len(records), "gzip_uncompressed_end": end_offset,
                  "counts": dict(current), "outputs": outputs, "completed_utc": stamp()}
            atomic_json(work / "qc.json", qc)
            for temporary_name in ("input.bed", "mapped.bed", "unmapped.bed"):
                (work / temporary_name).unlink(missing_ok=True)
            committed = chunks_root / f"chunk_{chunk_index:06d}"
            if committed.exists():
                raise RuntimeError(f"Unexpected concurrent checkpoint: {committed}")
            os.rename(work, committed)
            completed.append(committed)
            counts.update(current)
            offset += len(records)
            fraction = min(1.0, compressed.tell() / source_signature["size"])
            elapsed = time.monotonic() - started
            progress.update("chunk_committed", committed_rows=offset, gzip_percent=100 * fraction, counts=dict(counts))
            filled = int(fraction * 24)
            log(f"{trait} [{'#' * filled}{'-' * (24-filled)}] compressed~{100*fraction:5.1f}% "
                f"chunk={chunk_index} rows={offset:,} PASS={counts['PASS']:,} REVIEW={counts['REVIEW']:,} "
                f"chunk_seconds={elapsed:.1f} rate={len(records)/max(elapsed, 0.001):,.0f} rows/s")
            chunk_index += 1
            del records, destinations
    if not offset:
        raise RuntimeError(f"No data rows in {source}")
    if signature(source) != source_signature:
        raise RuntimeError(f"Input file changed while processing: {source}")
    if counts["INPUT"] != counts["PASS"] + counts["REVIEW"] or counts["INPUT"] != offset:
        raise RuntimeError("Final row conservation check failed")
    output_info = {}
    for name, key, header in ((f"{trait}.SNP.GRCh38.tsv.gz", "accepted", output_header),
                              (f"{trait}.review.tsv.gz", "review", output_header + ["raw_line_base64"])):
        output_info[name] = merged_file(directory / name, header, completed, key, progress)
    summary = {"version": VERSION, "config_id": config_id, "trait": trait, "counts": dict(counts),
               "outputs": output_info, "input_gzip_crc_checked_to_eof": True, "chunks": len(completed),
               "scope": "Validated biallelic SNPs; all other rows in review; target duplicates retained",
               "statistical_columns": "Original values unchanged; beta is not redefined relative to genomic ALT",
               "completed_utc": stamp()}
    atomic_json(done, summary)
    log(f"TRAIT_COMPLETE {trait}: INPUT={offset:,}, PASS={counts['PASS']:,}, REVIEW={counts['REVIEW']:,}")
    return summary


def qc_table(root, results):
    columns = ["trait", "INPUT", "PASS", "REVIEW", "ELIGIBLE_SNP", *STATUS[1:],
               "FASTA_CHECKED", "PLUS_CHECKED", "MINUS_CHECKED", "PALINDROMIC_PASS",
               "PAIR_MATCH_RATE", "PLUS_MATCH_RATE", "MINUS_MATCH_RATE"]
    temporary = root / "QC_summary.tsv.partial"
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for result in results:
            count = Counter(result["counts"])
            row = {key: count[key] for key in columns if key != "trait"}
            row["trait"] = result["trait"]
            for prefix in ("PAIR", "PLUS", "MINUS"):
                n = count["FASTA_CHECKED"] if prefix == "PAIR" else count[prefix + "_CHECKED"]
                passed = count["PASS"] if prefix == "PAIR" else count[prefix + "_PASS"]
                row[prefix + "_MATCH_RATE"] = f"{passed/n:.8f}" if n else "NA"
            writer.writerow(row)
    os.replace(temporary, root / "QC_summary.tsv")


def run(args):
    root = args.out
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".run.lock").open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another run already holds this output directory; no jobs were stopped")
        progress = Progress(root)
        fasta = None
        try:
            progress.update("preflight")
            paths = [args.lift, args.chain, args.fasta, Path(str(args.fasta) + ".fai")]
            sources = [args.base / f"data/GWAS/raw/{trait}_EA_GWAMA_20190214_GC.gz" for trait in TRAITS]
            for path in paths + sources:
                if not path.is_file():
                    raise RuntimeError(f"Required file not found: {path}")
            usage = subprocess.run([str(args.lift)], capture_output=True, text=True, timeout=15)
            if "usage:" not in (usage.stdout + usage.stderr).lower():
                raise RuntimeError(f"liftOver launch failed: {usage.stdout}{usage.stderr}")
            fasta = Fasta(args.fasta)
            source_sizes = inspect_chain(args.chain, fasta)
            config = {"version": VERSION, "chunk_size": args.chunk_size,
                      "source_files": [signature(path) for path in sources],
                      "chain": {**signature(args.chain), "sha256": sha256(args.chain)},
                      "liftOver": {**signature(args.lift), "sha256": sha256(args.lift)},
                      "fasta": signature(args.fasta), "fai": {**signature(paths[3]), "sha256": sha256(paths[3])}}
            config_id = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
            manifest = root / "run_manifest.json"
            if manifest.exists():
                if load_json(manifest)["config_id"] != config_id:
                    raise RuntimeError("Inputs/resources/chunk size differ from this run's manifest; use a new --out directory")
            else:
                existing = [p for p in root.iterdir() if p.name not in (".run.lock", "progress.json")]
                if existing:
                    raise RuntimeError("Output directory contains unrecognized files and has no run manifest")
                atomic_json(manifest, {"config_id": config_id, "configuration": config, "created_utc": stamp(),
                                       "script_sha256": sha256(Path(__file__)), "python": sys.version,
                                       "compression": args.pigz or "Python gzip"})
            log(f"VERSION={VERSION}; CHUNK={args.chunk_size:,}; pigz_threads={args.compress_threads}")
            log(f"OUTPUT={root}; free_disk={shutil.disk_usage(root).free/1024**3:,.1f} GiB")
            results = []
            for trait in TRAITS:
                result = process_trait(trait, args, root, config_id, source_sizes, fasta, progress)
                results.append(result)
                qc_table(root, results)
            atomic_json(root / "COMPLETE.json", {"config_id": config_id, "completed_utc": stamp(),
                        "traits": {result["trait"]: result["counts"] for result in results}})
            progress.update("COMPLETE", counts={r["trait"]: r["counts"] for r in results})
            log(f"ALL_THREE_COMPLETE | QC={root / 'QC_summary.tsv'}")
        except BaseException as error:
            progress.update("ERROR", error=f"{type(error).__name__}: {error}")
            raise
        finally:
            if fasta is not None:
                fasta.close()


def status(args):
    root = args.out
    print("OUTPUT:", root)
    if not root.is_dir():
        print("Not started: output directory is absent.")
        return
    if (root / "progress.json").is_file():
        print("\n=== Latest progress ===")
        print(json.dumps(load_json(root / "progress.json"), indent=2))
    print("\n=== Committed results ===")
    for trait in TRAITS:
        directory = root / trait
        counts, chunks = Counter(), 0
        for path in sorted((directory / "chunks").glob("chunk_[0-9]*/qc.json")):
            counts.update(load_json(path)["counts"])
            chunks += 1
        print(f"{trait}: chunks={chunks}, INPUT={counts['INPUT']:,}, PASS={counts['PASS']:,}, REVIEW={counts['REVIEW']:,}, "
              f"MINUS_CHECKED={counts['MINUS_CHECKED']:,}, NON_SNP={counts['NON_SNP_NEEDS_VARIANT_LIFTOVER']:,}, "
              f"ALLELE_MISMATCH={counts['ALLELE_MISMATCH']:,}, COMPLETE={(directory/'COMPLETE.json').is_file()}")
        for pattern in (f"{trait}.SNP.GRCh38.tsv.gz", f"{trait}.review.tsv.gz"):
            path = directory / pattern
            if path.is_file():
                print(f"  {path.name}: {path.stat().st_size/1024**3:.3f} GiB")
                if (directory / "COMPLETE.json").is_file():
                    with gzip.open(path, "rt") as handle:
                        for _ in range(2):
                            line = handle.readline().rstrip()
                            if line:
                                print("   ", line)
    summary = root / "QC_summary.tsv"
    if summary.is_file():
        print("\n=== QC summary ===")
        print(summary.read_text())
    print("ALL_THREE_COMPLETE:", (root / "COMPLETE.json").is_file())
    logs = list((args.base / "logs").glob("Step_03_GWAS_hg38_*.log"))
    if logs:
        latest = max(logs, key=lambda p: p.stat().st_mtime_ns)
        from collections import deque
        with latest.open(errors="replace") as handle:
            tail = deque(handle, maxlen=12)
        print("\n=== Latest log:", latest, "===")
        print("".join(tail))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "status"))
    parser.add_argument("--base", type=Path, default=Path.home() / "Drug_rankGWAS")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--atlas", type=Path, default=Path.home() / "GRN_GWAS_Atlas")
    parser.add_argument("--lift", type=Path)
    parser.add_argument("--chain", type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--chunk-size", type=int, default=500000)
    parser.add_argument("--compress-threads", type=int, default=2)
    parser.add_argument("--min-free-gib", type=float, default=5.0)
    args = parser.parse_args()
    args.base, args.atlas = args.base.expanduser().resolve(), args.atlas.expanduser().resolve()
    args.out = (args.out or args.base / "processed/GWAS/GRCh38_UCSC_v1").expanduser().resolve()
    args.lift = (args.lift or args.atlas / "tools/ucsc_liftover_env/bin/liftOver").expanduser().resolve()
    args.chain = (args.chain or args.base / "reference/liftover_ucsc/hg19ToHg38.over.chain.gz").expanduser().resolve()
    args.fasta = (args.fasta or args.atlas / "reference/GRCh38/GRCh38.primary_assembly.genome.fa").expanduser().resolve()
    args.pigz = shutil.which("pigz")
    if args.chunk_size < 1 or not 1 <= args.compress_threads <= 4 or args.min_free_gib < 0:
        parser.error("Positive chunk size, 1-4 compression threads, and a nonnegative disk reserve are required")
    try:
        (run if args.action == "run" else status)(args)
    except Exception:
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
