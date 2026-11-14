import contextlib
import gzip
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from workflow_config import parser, resolve


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


engine = module("production", "Step_03_liftover.py")
audit = module("final_audit", "Step_04_final_audit.py")
HEADERS = "rs_number reference_allele other_allele eaf beta se beta_95L beta_95U z p-value _-log10_p-value q_statistic q_p-value i2 n_studies n_samples effects chi-squared_GCcorr p-value_GCcorr".split()
COLUMNS = {key: i for i, key in enumerate(HEADERS)}


def record(marker="1:10_A_G", effect="G", other="A"):
    values = {key: "0.00123000" for key in HEADERS}
    values.update(rs_number=marker, reference_allele=effect, other_allele=other)
    return "\t".join(values[key] for key in HEADERS).encode()


class ConfigTests(unittest.TestCase):
    def test_profile_relative_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            profile = root / "profile.json"
            profile.write_text(json.dumps({"inputs":{"Test":"input.gz"},"out":"run"}))
            args = resolve(parser("test"), ["--config", str(profile)])
            self.assertEqual(args.inputs, {"Test": root / "input.gz"})
            self.assertEqual(args.out, root / "run")
            self.assertEqual(args.chunk_size, 500000)

    def test_cli_input_replaces_profile_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            profile = Path(tmp) / "profile.json"
            profile.write_text(json.dumps({"inputs":{"Old":"old.gz"}}))
            args = resolve(parser("test"), ["--config",str(profile),"--input","New=new.gz"])
            self.assertEqual(args.traits, ("New",))
            self.assertEqual(args.inputs["New"], Path.cwd() / "new.gz")

    def test_unsafe_trait_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            resolve(parser("test"), ["--input","../bad=anything.gz"])

    def test_duplicate_trait_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            resolve(parser("test"), ["--input","A=1.gz","--input","A=2.gz"])

    def test_unknown_profile_key_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            profile = Path(tmp) / "profile.json"
            profile.write_text('{"chunk_szie": 100}')
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                resolve(parser("test"),["--config",str(profile)])


class ParsingTests(unittest.TestCase):
    def parse(self, raw):
        return engine.parse_record(raw, 7, COLUMNS, b"\t", {"chr1":1000})

    def test_supported_snp(self):
        result = self.parse(record())
        self.assertEqual(result[1:4], (7,"chr1",10))
        self.assertEqual(result[-1], "NO_MAP")
        self.assertEqual(result[0], record())

    def test_marker_order_is_not_effect_order(self):
        self.assertEqual(self.parse(record("1:10_G_A","A","G"))[-1], "NO_MAP")

    def test_indel_is_reviewed(self):
        self.assertEqual(self.parse(record("1:10_A_AT","AT","A"))[-1], "NON_SNP_NEEDS_VARIANT_LIFTOVER")

    def test_allele_set_mismatch(self):
        self.assertEqual(self.parse(record("1:10_A_G","A","C"))[-1], "ID_ALLELE_MISMATCH")

    def test_bad_row_and_id(self):
        self.assertEqual(self.parse(b"bad\trow")[-1], "BAD_COLUMN_COUNT")
        self.assertEqual(self.parse(record("rs123"))[-1], "BAD_ID")

    def test_out_of_bounds(self):
        self.assertEqual(self.parse(record("1:1001_A_G"))[-1], "BAD_SOURCE_COORD")


class OutputTests(unittest.TestCase):
    def test_reverse_complement_keeps_statistical_text(self):
        class Fasta:
            def base(self, chrom, pos):
                return b"T"
        parsed = engine.parse_record(record("1:10_A_C","C","A"),1,COLUMNS,b"\t",{"chr1":1000})
        class Progress:
            def update(self, *a, **kw):
                pass
        with tempfile.TemporaryDirectory() as tmp:
            counts, accepted, review = engine.write_rows([parsed],{1:("chr1",20,21,"-")},Fasta(),Path(tmp),b"\t",19,Progress())
            fields = accepted.read_bytes().rstrip(b"\n").split(b"\t")
            self.assertEqual(fields[:19], record("1:10_A_C","C","A").split(b"\t"))
            added = dict(zip(engine.EXTRA,fields[19:]))
            self.assertEqual(added["reference_allele_hg38"], b"G")
            self.assertEqual(added["other_allele_hg38"], b"T")
            self.assertEqual(added["REF_hg38"],b"T")
            self.assertEqual(added["ALT_hg38"],b"G")
            self.assertEqual(added["ALT_hg38_from"],b"reference_allele_GRCh37")
            self.assertEqual(counts["PASS"],1)
            self.assertEqual(counts["MINUS_PASS"],1)
            self.assertEqual(review.read_bytes(),b"")

    def test_final_file_hash_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"part.gz"
            path.write_bytes(gzip.compress(b"one\ntwo\n"))
            info = engine.file_info(path)
            engine.verify(path,info)
            data = bytearray(path.read_bytes());data[-1] ^= 1;path.write_bytes(data)
            with self.assertRaises(RuntimeError):
                engine.verify(path,info)

    def test_row_accounting_failure(self):
        c = {"INPUT":1,"PASS":1,"REVIEW":0,"ELIGIBLE_SNP":1,"FASTA_CHECKED":1,
             "PLUS_CHECKED":1,"PLUS_PASS":1}
        audit.check_counts(c,"test")
        c["INPUT"] = 2
        with self.assertRaises(RuntimeError):
            audit.check_counts(c,"test")

    def test_recorded_production_counts_are_internally_consistent(self):
        data = json.loads((ROOT/"results/production_qc.json").read_text())
        for item in data["traits"]:
            audit.check_counts(item["counts"],item["trait"])
            audit.check_summary_row(item["summary"],Counter(item["counts"]),item["trait"])


if __name__ == "__main__":
    unittest.main()
