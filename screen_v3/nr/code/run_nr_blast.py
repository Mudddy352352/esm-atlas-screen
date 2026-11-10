#!/usr/bin/env python3
"""The real nr search, which section 3.8 has so far only projected.

Section 3.8 reports a BLAST-class aligner run on local sequence sets with its
E-values scaled to the size of a large database, because no sequence search
service was reachable when the screen was run. NCBI's BLAST URL API answers
from this environment now, so this runs the search itself: each shortlisted
Atlas protein against nr, and the question the projection was standing in for,
which is whether a pairwise search links the protein to a characterised
drug-metabolising enzyme.

One job at a time, polled politely. Resumable: request ids and finished
results are written as they arrive, so a restart picks up where it stopped and
never resubmits a job that is already running.

Usage: run_nr_blast.py <work_dir> [max_minutes]
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
TIER1 = _ROOT + "/screen_v2/results/tier1_candidates.tsv"
FASTA = _ROOT + "/screen_v2/taxonomy/candidates.fasta"
TOOL = "esm-atlas-screen"
POLL = 60


def post(params, timeout=180):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(BLAST, data=data)
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf8", "replace")


def sequences():
    out, name, buf = {}, None, []
    for line in open(FASTA):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def shortlist():
    rows = [l.rstrip("\n").split("\t") for l in open(TIER1)]
    head = rows[0]
    return [dict(zip(head, r)) for r in rows[1:]]


def submit(seq):
    txt = post({"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
                "QUERY": seq, "HITLIST_SIZE": "100", "EXPECT": "10",
                "TOOL": TOOL})
    rid = re.search(r"RID = (\S+)", txt)
    if not rid:
        raise RuntimeError("no RID returned")
    return rid.group(1)


def ready(rid):
    txt = post({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid,
                "TOOL": TOOL}, timeout=120)
    if "Status=WAITING" in txt:
        return "waiting"
    if "Status=FAILED" in txt or "Status=UNKNOWN" in txt:
        return "failed"
    if "Status=READY" in txt:
        return "no hits" if "ThereAreHits=yes" not in txt else "ready"
    return "waiting"


def fetch(rid):
    return post({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
                 "ALIGNMENTS": "100", "DESCRIPTIONS": "100",
                 "TOOL": TOOL}, timeout=300)


def main(work, max_minutes=600):
    max_minutes = float(max_minutes)
    res = os.path.join(work, "results")
    os.makedirs(res, exist_ok=True)
    state_path = os.path.join(res, "state.json")
    state = json.load(open(state_path)) if os.path.exists(state_path) else {}
    seqs = sequences()
    started = time.time()

    for row in shortlist():
        acc = row["target"]
        if state.get(acc, {}).get("done"):
            continue
        if (time.time() - started) / 60 > max_minutes:
            print("time limit reached, %d left" % sum(
                1 for r in shortlist() if not state.get(r["target"], {}).get("done")),
                flush=True)
            break
        seq = seqs.get(acc)
        if not seq:
            state[acc] = dict(done=True, note="no sequence in the candidate fasta")
            json.dump(state, open(state_path, "w"), indent=1)
            continue
        rid = state.get(acc, {}).get("rid")
        if not rid:
            rid = submit(seq)
            state[acc] = dict(rid=rid, done=False, query=row["query"])
            json.dump(state, open(state_path, "w"), indent=1)
            print("[%s] submitted, RID %s" % (acc, rid), flush=True)
            time.sleep(10)
        while True:
            st = ready(rid)
            if st != "waiting":
                break
            if (time.time() - started) / 60 > max_minutes:
                print("[%s] still running at the time limit, RID kept" % acc, flush=True)
                return
            time.sleep(POLL)
        if st == "ready":
            open(os.path.join(res, "%s.xml" % acc), "w").write(fetch(rid))
        state[acc] = dict(rid=rid, done=True, status=st, query=row["query"])
        json.dump(state, open(state_path, "w"), indent=1)
        print("[%s] %s" % (acc, st), flush=True)
        time.sleep(5)
    print("finished pass", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
