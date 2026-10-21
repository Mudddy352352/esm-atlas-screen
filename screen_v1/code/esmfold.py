#!/usr/bin/env python3
"""Fold a FASTA sequence with the public ESMFold API (api.esmatlas.com).

The API accepts single sequences up to 400 residues and is rate/time limited,
so each request is retried. Writes a PDB next to the output directory.

Usage: esmfold.py <outdir> <query.fasta> [more.fasta ...]
       esmfold.py --window START END <outdir> <query.fasta>   (1-based, inclusive)
"""
import os
import sys
import time
import urllib.error
import urllib.request

URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAXLEN = 400


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            seq.append(line.strip())
    return name, "".join(seq)


def fold(seq, tries=6):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(URL, data=seq.encode(),
                                         headers={"Content-Type": "text/plain"})
            with urllib.request.urlopen(req, timeout=600) as r:
                body = r.read().decode()
            if body.lstrip().startswith(("HEADER", "ATOM", "PARENT")):
                return body
            last = body[:200]
        except Exception as e:
            last = repr(e)
        wait = 20 * (i + 1)
        print(f"    retry {i+1}/{tries} after {wait}s ({last})", flush=True)
        time.sleep(wait)
    raise RuntimeError(f"ESMFold failed: {last}")


def main():
    args = sys.argv[1:]
    window = None
    if args[0] == "--window":
        window = (int(args[1]), int(args[2]))
        args = args[3:]
    outdir, paths = args[0], args[1:]
    os.makedirs(outdir, exist_ok=True)
    for path in paths:
        base = os.path.basename(path).replace(".fasta", "")
        name, seq = read_fasta(path)
        if window:
            seq = seq[window[0] - 1:window[1]]
            base += f"_w{window[0]}-{window[1]}"
        out = os.path.join(outdir, base + ".pdb")
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            print(f"[skip] {base}", flush=True)
            continue
        if len(seq) > MAXLEN:
            print(f"[skip] {base}: {len(seq)} aa exceeds API limit {MAXLEN}", flush=True)
            continue
        print(f"[fold] {base} ({len(seq)} aa)", flush=True)
        pdb = fold(seq)
        with open(out, "w") as fh:
            fh.write(pdb)
        # mean pLDDT from the B-factor column of CA atoms
        bs = [float(l[60:66]) for l in pdb.splitlines()
              if l.startswith("ATOM") and l[12:16].strip() == "CA"]
        print(f"[done] {base} residues={len(bs)} mean_pLDDT={sum(bs)/len(bs):.1f}", flush=True)
        time.sleep(5)


if __name__ == "__main__":
    main()
