#!/usr/bin/env python3
"""Fetch the official forward chain and optionally the Linux x86_64 liftOver binary.

External resources are downloaded at run time and are not bundled in releases.
Their upstream terms apply: https://genome.ucsc.edu/license/
"""
import argparse
import gzip
import hashlib
import os
import tempfile
import urllib.request
from pathlib import Path


def fetch(url, destination, binary=False):
    if destination.exists():
        raise FileExistsError(f"Preserving existing resource: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=".download_", dir=destination.parent)
    os.close(descriptor)
    temporary = Path(name)
    try:
        print("[DOWNLOAD]", url, flush=True)
        with urllib.request.urlopen(url, timeout=30) as source, temporary.open("wb") as out:
            total = 0
            while True:
                block = source.read(1024**2)
                if not block:
                    break
                out.write(block)
                total += len(block)
                if total > 128 * 1024**2:
                    raise RuntimeError("Unexpected resource size")
        if binary:
            with temporary.open("rb") as handle:
                if handle.read(4) != b"\x7fELF":
                    raise RuntimeError("Expected a Linux ELF executable")
            temporary.chmod(0o755)
        else:
            with gzip.open(temporary, "rt", encoding="ascii") as handle:
                headers = [line.split() for line in handle if line.startswith("chain ")]
            if not any(len(h)==13 and h[2:5]==["chr1","249250621","+"] and h[7:9]==["chr1","248956422"] for h in headers):
                raise RuntimeError("Forward hg19-to-hg38 chain check failed")
        os.replace(temporary, destination)
        print("[SAVED]", destination, "SHA256", hashlib.sha256(destination.read_bytes()).hexdigest(), flush=True)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chain", type=Path, default=Path("reference/liftover/hg19ToHg38.over.chain.gz"))
    p.add_argument("--binary", type=Path, help="Also fetch Linux x86_64 liftOver to this path")
    args = p.parse_args()
    if not args.chain.exists():
        fetch("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", args.chain)
    else:
        print("[EXISTS]", args.chain)
    if args.binary:
        if args.binary.exists():
            print("[EXISTS]", args.binary)
        else:
            fetch("https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver", args.binary, binary=True)


if __name__ == "__main__":
    main()
