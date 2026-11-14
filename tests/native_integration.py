import contextlib
import csv
import fcntl
import gzip
import importlib.util
import io
import json
import os
import sys
import subprocess
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))
spec = importlib.util.spec_from_file_location("pipeline", HERE / "scripts/Step_03_liftover.py")
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)
LIFT = Path(os.environ.get("LIFTOVER_BIN", str(HERE / "tools/liftOver"))).resolve()
if not LIFT.is_file():
    raise SystemExit("Set LIFTOVER_BIN to a local UCSC liftOver executable; none is bundled.")
NAMES = "rs_number reference_allele other_allele eaf beta se beta_95L beta_95U z p-value _-log10_p-value q_statistic q_p-value i2 n_studies n_samples effects chi-squared_GCcorr p-value_GCcorr".split()


def write_fasta(path):
    one = list("A" * 200)
    one[20], one[22], one[169], one[80], one[90] = "A", "A", "T", "C", "N"
    index = []
    with path.open("wb") as handle:
        for name, seq in (("chr1", "".join(one)), ("chr2", "C" * 200), ("chrUn", "A" * 10)):
            handle.write((">" + name + "\n").encode())
            offset = handle.tell()
            for i in range(0, len(seq), 7):
                handle.write((seq[i:i+7].lower() + "\n").encode())
            index.append(f"{name}\t{len(seq)}\t{offset}\t7\t8\n")
    Path(str(path) + ".fai").write_text("".join(index))


def chain_text():
    entries = [(0, 10, "chr1", 200, "+", 20, 30),
               (20, 30, "chr1", 200, "-", 30, 40),
               (40, 41, "chr1", 200, "+", 50, 51),
               (40, 41, "chr2", 200, "+", 60, 61),
               (60, 61, "chr1", 200, "+", 80, 81),
               (70, 71, "chr1", 200, "+", 90, 91),
               (80, 81, "chrUn", 10, "+", 0, 1)]
    return "\n".join(f"chain 1000 chr1 1000 + {s} {e} {q} {size} {strand} {qs} {qe} {i}\n{e-s}\n"
                     for i, (s,e,q,size,strand,qs,qe) in enumerate(entries,1)) + "\n"


def rows():
    variants = [("1:1_A_G", "G", "A"), ("1:21_A_C", "C", "A"),
                ("1:41_A_G", "G", "A"), ("1:51_A_G", "G", "A"),
                ("1:55_A_AT", "AT", "A"), ("1:61_A_G", "G", "A"),
                ("1:71_A_G", "G", "A"), ("1:81_A_G", "G", "A"),
                ("broken", "G", "A"), ("1:2_A_G", "C", "T")]
    result = []
    for marker, ra, oa in variants:
        row = {key: "0.001230000" for key in NAMES}
        row.update(rs_number=marker, reference_allele=ra, other_allele=oa,
                   beta="-0.0012300", eaf="0.43353100000000006", effects="--++-+", n_samples="486767")
        result.append("\t".join(row[k] for k in NAMES))
    result.append("malformed\tline")
    row = {key: "1e-20" for key in NAMES}
    row.update(rs_number="1:3_A_T", reference_allele="T", other_allele="A", beta="-0.0012300", eaf="0.43353100000000006")
    result.append("\t".join(row[k] for k in NAMES))
    return result


def execute(args):
    with contextlib.redirect_stdout(io.StringIO()):
        p.run(args)


with tempfile.TemporaryDirectory(prefix="gwas_native_test_") as temporary:
    root = Path(temporary)
    base = root / "Drug_rankGWAS"
    raw = base / "data/GWAS/raw"
    raw.mkdir(parents=True)
    fasta = root / "reference.fa"
    write_fasta(fasta)
    chain = root / "toy.chain.gz"
    with gzip.open(chain, "wt") as handle:
        handle.write(chain_text())
    for trait in p.TRAITS:
        with gzip.open(raw / f"{trait}_EA_GWAMA_20190214_GC.gz", "wt") as handle:
            handle.write("\t".join(NAMES) + "\n" + "\n".join(rows()) + "\n")
    args = SimpleNamespace(base=base, out=root / "output", lift=LIFT, chain=chain, fasta=fasta,
                           chunk_size=4, pigz=shutil.which("pigz"), compress_threads=2, min_free_gib=0.05,
                           traits=p.TRAITS, inputs={t: raw / f"{t}_EA_GWAMA_20190214_GC.gz" for t in p.TRAITS})
    # Use miniature genomes for native liftOver tests, bypassing only the human-assembly preflight.
    original_inspect = p.inspect_chain
    p.inspect_chain = lambda chain, fasta: {"chr1": 1000}
    original_map = p.map_chunk
    calls = []
    def interrupt_after_one(records, *a):
        calls.append(records[0][1])
        if len(calls) == 2:
            raise RuntimeError("TEST interruption after first committed chunk")
        return original_map(records, *a)
    p.map_chunk = interrupt_after_one
    try:
        execute(args)
        raise AssertionError("Expected interruption")
    except RuntimeError as error:
        assert "TEST interruption" in str(error)
    first = args.out / "MCH/chunks/chunk_000000"
    original_stat = (first / "accepted.tsv.gz").stat().st_mtime_ns
    assert first.is_dir() and not (args.out / "COMPLETE.json").exists()
    calls = []
    def track_resume(records, *a):
        calls.append(records[0][1])
        return original_map(records, *a)
    p.map_chunk = track_resume
    execute(args)
    assert calls[:2] == [5,9], calls
    assert (first / "accepted.tsv.gz").stat().st_mtime_ns == original_stat
    assert (args.out / "COMPLETE.json").is_file()
    expected = {"INPUT":12,"PASS":3,"REVIEW":9,"ELIGIBLE_SNP":8,"NO_MAP":1,"MULTI_MAP":1,
                "NON_SNP_NEEDS_VARIANT_LIFTOVER":1,"ALLELE_MISMATCH":1,"FASTA_N":1,
                "NON_PRIMARY_OR_NON_1BP":1,"BAD_ID":1,"ID_ALLELE_MISMATCH":1,
                "BAD_COLUMN_COUNT":1,"PLUS_PASS":2,"MINUS_PASS":1,"PALINDROMIC_PASS":1}
    for trait in p.TRAITS:
        directory = args.out / trait
        summary = json.loads((directory / "COMPLETE.json").read_text())
        for key, value in expected.items():
            assert summary["counts"][key] == value, (trait,key,summary["counts"])
        with gzip.open(directory / f"{trait}.SNP.GRCh38.tsv.gz", "rt") as handle:
            accepted = list(csv.DictReader(handle, delimiter="\t"))
        with gzip.open(directory / f"{trait}.review.tsv.gz", "rt") as handle:
            review = list(csv.DictReader(handle, delimiter="\t"))
        assert len(accepted) == 3 and len(review) == 9
        assert accepted[0]["pos38"] == "21" and accepted[0]["REF_hg38"] == "A"
        assert accepted[1]["pos38"] == "170" and accepted[1]["lift_strand"] == "-"
        assert (accepted[1]["reference_allele_hg38"],accepted[1]["other_allele_hg38"]) == ("G","T")
        assert accepted[1]["REF_hg38"] == "T" and accepted[1]["ALT_hg38"] == "G"
        assert accepted[1]["ALT_hg38_from"] == "reference_allele_GRCh37"
        assert all(row["beta"] == "-0.0012300" and row["eaf"] == "0.43353100000000006" for row in accepted)
        assert any(row["liftover_qc"] == "BAD_COLUMN_COUNT" and row["raw_line_base64"] != "NA" for row in review)
        assert {int(row["source_row"]) for row in accepted+review} == set(range(1,13))
    # Completed output reuse must not run any mapping again.
    p.map_chunk = lambda *a: (_ for _ in ()).throw(AssertionError("Unexpected recalculation"))
    execute(args)
    p.map_chunk = original_map
    # A lock held by another process/file description must prevent a second run.
    with (args.out / ".run.lock").open("a+") as held:
        fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            execute(args)
            raise AssertionError("Expected locking failure")
        except RuntimeError as error:
            assert "already holds" in str(error)
    # Corrupt a committed checkpoint and require safe refusal on resume.
    args.out = root / "corrupt_checkpoint"
    calls = []
    p.map_chunk = interrupt_after_one
    try:
        execute(args)
    except RuntimeError:
        pass
    corrupted = args.out / "MCH/chunks/chunk_000000/accepted.tsv.gz"
    data = bytearray(corrupted.read_bytes())
    data[-1] ^= 1
    corrupted.write_bytes(data)
    p.map_chunk = original_map
    try:
        execute(args)
        raise AssertionError("Expected checkpoint integrity failure")
    except RuntimeError as error:
        assert "integrity" in str(error)
    assert corrupted.read_bytes() == data
    # Source CRC failures must prevent completion, including on the last gzip member.
    p.map_chunk = original_map
    args.out = root / "source_crc_failure"
    source = raw / "MCH_EA_GWAMA_20190214_GC.gz"
    original_bytes = source.read_bytes()
    original_stat = source.stat()
    bad_source = bytearray(original_bytes)
    bad_source[-8] ^= 1
    source.write_bytes(bad_source)
    try:
        execute(args)
        raise AssertionError("Expected source CRC failure")
    except gzip.BadGzipFile:
        pass
    assert not (args.out / "COMPLETE.json").exists()
    source.write_bytes(original_bytes)
    os.utime(source, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    # Resource/input identity guards must reject mixed run configurations.
    args.out = root / "output"
    args.chunk_size = 5
    try:
        execute(args)
        raise AssertionError("Expected configuration guard failure")
    except RuntimeError as error:
        assert "manifest" in str(error)
    args.chunk_size = 4
    # Check a custom trait label and the portable final auditor on actual native outputs.
    args.out = root / "custom_trait"
    args.traits = ("Example",)
    args.inputs = {"Example": raw / "MCH_EA_GWAMA_20190214_GC.gz"}
    execute(args)
    audit_spec = importlib.util.spec_from_file_location("public_audit", HERE / "scripts/Step_04_final_audit.py")
    auditor = importlib.util.module_from_spec(audit_spec)
    audit_spec.loader.exec_module(auditor)
    with contextlib.redirect_stdout(io.StringIO()):
        report = auditor.audit(args.out)
    assert report["audit"] == "SNP_LIFTOVER_FINAL_AUDIT_PASS"
    assert report["results"][0]["trait"] == "Example"
    # Exercise the real assembly guard with chain headers and an independent size checker.
    p.inspect_chain = original_inspect
    human_sizes = {c: (249250621,248956422) if c=="chr1" else (135534747,133797422) if c=="chr10" else (1000,900) for c in p.PRIMARY}
    headers = root / "headers.chain.gz"
    with gzip.open(headers,"wt") as handle:
        for i,(chrom,(old,new)) in enumerate(human_sizes.items(),1):
            handle.write(f"chain 1 {chrom} {old} + 0 1 {chrom} {new} + 0 1 {i}\n1\n\n")
    class SizeChecker:
        def check_length(self, chrom, length):
            assert length == human_sizes[chrom][1]
    original_inspect(headers, SizeChecker())
    reader = p.Fasta(fasta)
    assert reader.base("chr1",7) == b"A" and reader.base("chr1",8) == b"A"
    assert reader.base("chr1",170) == b"T"
    reader.close()
    # Exercise the format-specific prefix pilot and inspection CLI using native mappings.
    pilot_spec = importlib.util.spec_from_file_location("public_pilot", HERE / "scripts/Step_02_pilot.py")
    pilot_module = importlib.util.module_from_spec(pilot_spec)
    pilot_spec.loader.exec_module(pilot_module)
    pilot_output = root / "pilot_test"
    pilot_output.mkdir()
    reader = p.Fasta(fasta)
    class PilotFasta:
        def base(self, chrom, position):
            return reader.base(chrom, position).decode("ascii")
    with contextlib.redirect_stdout(io.StringIO()):
        prefix = pilot_module.pilot("Example", args.inputs["Example"], 10, LIFT, chain, PilotFasta(), pilot_output)
    reader.close()
    assert prefix["SAMPLED"] == 10 and prefix["PASS"] == 2 and prefix["MINUS_CHECKED"] == 1
    inspected = subprocess.run([sys.executable, str(HERE / "scripts/Step_01_inspect.py"),
        "--input", "Example=" + str(args.inputs["Example"]), "--liftover", str(LIFT),
        "--chain", str(chain), "--fasta", str(fasta)], check=True, capture_output=True, text=True)
    assert json.loads(inspected.stdout)["inputs"]["Example"]["column_count"] == 19

print("PASS: actual UCSC liftOver on plus/minus synthetic chains; multiple/no mappings; target-base mismatches; indel and malformed-row audit; exact beta/eaf preservation; row conservation; concatenated gzip outputs; interrupted-run resume without recalculation; completed-output reuse; checkpoint corruption rejection; source CRC and configuration guards; duplicate-run locking; assembly direction and FASTA offset guards; inspection CLI, prefix pilot and custom-trait final audit.")
print("The user's HPC data were not accessed. Human assembly preflight was independently tested; miniature-genome integration bypassed only that preflight.")
