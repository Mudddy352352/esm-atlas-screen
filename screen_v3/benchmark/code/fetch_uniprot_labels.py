#!/usr/bin/env python3
"""Pull EC numbers and review status for the AlphaFold DB targets of the
labelled-hit benchmark, so that ground truth comes from UniProt curation
rather than from the free-text name in the Foldseek output.

Records are cached under labels/ so the step reproduces without network.

Usage: fetch_uniprot_labels.py <scored_tsv> <cache_dir> <out_tsv>
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

API = "https://rest.uniprot.org"
FIELDS = "accession,reviewed,protein_name,ec,xref_pfam,organism_name"
CHUNK = 2000


def post(url, data):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode(),
                                 method="POST")
    with urllib.request.urlopen(req, timeout=120) as fh:
        return json.load(fh)


def get(url):
    with urllib.request.urlopen(url, timeout=300) as fh:
        return fh.read().decode()


def run_chunk(accs, cache_path):
    if os.path.exists(cache_path):
        return open(cache_path).read()
    job = post(API + "/idmapping/run",
               {"from": "UniProtKB_AC-ID", "to": "UniProtKB", "ids": ",".join(accs)})
    jid = job["jobId"]
    for _ in range(120):
        st = json.loads(get(API + "/idmapping/status/" + jid))
        if st.get("jobStatus") in (None, "FINISHED") and ("results" in st or st.get("jobStatus") == "FINISHED"):
            break
        if st.get("jobStatus") == "ERROR":
            raise RuntimeError(st)
        time.sleep(3)
    txt = get("%s/idmapping/uniprotkb/results/stream/%s?fields=%s&format=tsv"
              % (API, jid, FIELDS))
    with open(cache_path, "w") as fh:
        fh.write(txt)
    return txt


def main(scored_tsv, cache_dir, out_tsv):
    os.makedirs(cache_dir, exist_ok=True)
    d = pd.read_csv(scored_tsv, sep="\t")
    afdb = d[d.db == "afdb50"].target.dropna().unique().tolist()
    accs = sorted({m.group(1) for t in afdb
                   for m in [re.match(r"AF-([A-Z0-9]+)-F\d+", t)] if m})
    print("unique AlphaFold accessions:", len(accs), flush=True)

    frames = []
    for i in range(0, len(accs), CHUNK):
        chunk = accs[i:i + CHUNK]
        path = os.path.join(cache_dir, "uniprot_%05d.tsv" % i)
        txt = run_chunk(chunk, path)
        frames.append(pd.read_csv(pd.io.common.StringIO(txt), sep="\t"))
        print("chunk %d: %d ids, %d rows" % (i, len(chunk), len(frames[-1])), flush=True)

    lab = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["Entry"])
    lab = lab.rename(columns={"Entry": "uniprot", "Reviewed": "reviewed",
                              "Protein names": "protein_name", "EC number": "ec",
                              "Pfam": "pfam", "Organism": "organism"})
    lab.to_csv(out_tsv, sep="\t", index=False)
    print("wrote", out_tsv, lab.shape, flush=True)
    print("with an EC:", int(lab.ec.notna().sum()), " reviewed:",
          int((lab.reviewed == "reviewed").sum()), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
