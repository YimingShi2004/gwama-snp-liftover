"""Shared path and input configuration for the portable command-line scripts."""
import argparse
import json
import re
import shutil
from pathlib import Path

DEFAULT_INPUTS = {trait: f"data/GWAS/raw/{trait}_EA_GWAMA_20190214_GC.gz"
                  for trait in ("MCH", "MCV", "RDW")}


def parser(description, action=False, pilot=False):
    p = argparse.ArgumentParser(description=description)
    if action:
        p.add_argument("action", choices=("run", "status"))
    p.add_argument("--config", type=Path, help="JSON profile; its paths are relative to the profile directory")
    p.add_argument("--input", action="append", metavar="TRAIT=FILE", help="Repeat to replace the profile's input mapping")
    p.add_argument("--out", type=Path, help="Output directory; command-line paths are relative to the current directory")
    p.add_argument("--liftover", "--lift", dest="lift", type=Path)
    p.add_argument("--chain", type=Path)
    p.add_argument("--fasta", type=Path)
    p.add_argument("--chunk-size", type=int)
    p.add_argument("--compress-threads", type=int)
    p.add_argument("--min-free-gib", type=float)
    if pilot:
        p.add_argument("--rows", type=int)
    return p


def resolve(p, argv=None):
    args = p.parse_args(argv)
    config, base = {}, Path.cwd()
    if args.config:
        profile = args.config.expanduser().resolve()
        base = profile.parent
        with profile.open() as handle:
            config = json.load(handle)
        if not isinstance(config, dict):
            p.error("The profile must be a JSON object")
        unknown = set(config) - {"inputs", "out", "liftover", "chain", "fasta", "chunk_size",
                                 "compress_threads", "min_free_gib", "pilot_rows"}
        if unknown:
            p.error("Unknown profile keys: " + ", ".join(sorted(unknown)))

    def path(value, relative_to):
        value = Path(value).expanduser()
        return (value if value.is_absolute() else relative_to / value).resolve()

    inputs = config.get("inputs", DEFAULT_INPUTS)
    if args.input:
        inputs = {}
        for item in args.input:
            if "=" not in item:
                p.error("--input must have the form TRAIT=FILE")
            trait, filename = item.split("=", 1)
            if trait in inputs:
                p.error("Duplicate input trait: " + trait)
            inputs[trait] = path(filename, Path.cwd())
    if not isinstance(inputs, dict) or not inputs:
        p.error("At least one named input is required")
    for trait, filename in inputs.items():
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,63}", trait) or not str(filename):
            p.error("Input traits must be simple labels of at most 64 characters")
    args.inputs = {trait: path(filename, base) for trait, filename in inputs.items()}
    args.traits = tuple(args.inputs)
    defaults = {"out": "results/run", "lift": "tools/liftOver",
                "chain": "reference/liftover/hg19ToHg38.over.chain.gz",
                "fasta": "reference/GRCh38/GRCh38.primary_assembly.genome.fa"}
    for key, default in defaults.items():
        given = getattr(args, key)
        config_key = "liftover" if key == "lift" else key
        setattr(args, key, path(given, Path.cwd()) if given is not None else path(config.get(config_key, default), base))
    for key, default in (("chunk_size", 500000), ("compress_threads", 2), ("min_free_gib", 5.0)):
        given = getattr(args, key)
        setattr(args, key, given if given is not None else config.get(key, default))
    if not isinstance(args.chunk_size, int) or args.chunk_size < 1:
        p.error("chunk_size must be a positive integer")
    if not isinstance(args.compress_threads, int) or not 1 <= args.compress_threads <= 4:
        p.error("compress_threads must be an integer from 1 to 4")
    if not isinstance(args.min_free_gib, (int, float)) or args.min_free_gib < 0:
        p.error("min_free_gib must be nonnegative")
    if hasattr(args, "rows"):
        args.rows = args.rows if args.rows is not None else config.get("pilot_rows", 10000)
        if not isinstance(args.rows, int) or args.rows < 1:
            p.error("pilot_rows must be a positive integer")
    args.base, args.pigz = base, shutil.which("pigz")
    return args
