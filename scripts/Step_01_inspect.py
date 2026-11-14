#!/usr/bin/env python3
"""Preview named gzip GWAMA inputs and the configured local resources."""
import gzip
import json
import shutil
from workflow_config import parser, resolve


def main(argv=None):
    args = resolve(parser(__doc__), argv)
    preview = {}
    for trait, source in args.inputs.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        with gzip.open(source, "rt", encoding="utf-8-sig") as handle:
            header = handle.readline().rstrip("\r\n")
            separator = "\t" if "\t" in header else None
            names = header.split(separator)
            rows = [handle.readline().rstrip("\r\n") for _ in range(2)]
        preview[trait] = {"file": str(source), "gzip_bytes": source.stat().st_size, "columns": names,
                          "column_count": len(names), "example_rows": rows}
    resources = {key: {"path": str(getattr(args, key)), "exists": getattr(args, key).is_file()}
                 for key in ("lift", "chain", "fasta")}
    resources["fai"] = {"path": str(args.fasta) + ".fai", "exists": args.fasta.with_name(args.fasta.name + ".fai").is_file()}
    print(json.dumps({"inputs": preview, "resources": resources, "pigz": shutil.which("pigz"),
                      "scope": "Prefix inspection only; does not establish the source genome assembly or full gzip integrity"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
