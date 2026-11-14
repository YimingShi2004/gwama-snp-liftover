#!/usr/bin/env python3
"""Export run provenance with absolute paths replaced by basenames."""
import argparse
import json
from pathlib import Path


def redact(value):
    if isinstance(value, dict):
        return {key: redact(item) for key, item in value.items()}
    if isinstance(value, list):
        return [redact(item) for item in value]
    if isinstance(value, str) and value.startswith("/"):
        return Path(value).name
    return value


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    with (args.root / "run_manifest.json").open() as handle:
        manifest = json.load(handle)
    with (args.root / "COMPLETE.json").open() as handle:
        complete = json.load(handle)
    traits = manifest["configuration"].get("traits", list(complete["traits"]))
    results = {}
    for trait in traits:
        with (args.root / trait / "COMPLETE.json").open() as handle:
            results[trait] = json.load(handle)
    data = redact({"run_manifest": manifest, "completion": complete, "traits": results})
    data["export_note"] = "Absolute path strings were replaced by basenames; retained config_id identifies the original unredacted configuration. Review before publishing."
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(data, indent=2) + "\n")
    print(args.out)


if __name__ == "__main__":
    main()
