#!/usr/bin/env python3
"""Submit every outstanding shortlisted protein to NCBI nr and record its id.

Background processes do not survive between turns in this environment, so the
search cannot be driven by a long-running poll loop. NCBI holds a finished
result for 24 hours under its request id, so the work is split: this submits
everything and writes the ids down, and collect_nr.py fetches whatever is ready
whenever it is next run. Nothing is resubmitted if it already has an id.

Usage: submit_all.py <work_dir>
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_nr_blast import BLAST, TOOL, sequences, shortlist  # noqa: E402


def submit(seq):
    d = {"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr", "QUERY": seq,
         "HITLIST_SIZE": "100", "EXPECT": "10", "TOOL": TOOL}
    req = urllib.request.Request(BLAST, data=urllib.parse.urlencode(d).encode())
    txt = urllib.request.urlopen(req, timeout=180).read().decode("utf8", "replace")
    m = re.search(r"RID = (\S+)", txt)
    if not m:
        raise RuntimeError("no RID returned")
    return m.group(1)


def main(work):
    res = os.path.join(work, "results")
    path = os.path.join(res, "state.json")
    state = json.load(open(path)) if os.path.exists(path) else {}
    seqs = sequences()
    for row in shortlist():
        acc = row["target"]
        cur = state.get(acc, {})
        if cur.get("done") or cur.get("rid"):
            continue
        seq = seqs.get(acc)
        if not seq:
            state[acc] = dict(done=True, note="no sequence in the candidate fasta")
            continue
        try:
            rid = submit(seq)
            state[acc] = dict(rid=rid, done=False, query=row["query"],
                              submitted=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            print("[%s] %s" % (acc, rid), flush=True)
        except Exception as exc:
            print("[%s] submit failed: %s" % (acc, exc), flush=True)
        json.dump(state, open(path, "w"), indent=1)
        time.sleep(12)
    n = sum(1 for v in state.values() if v.get("done"))
    print("ids held: %d, finished: %d, of %d" % (len(state), n, len(shortlist())))


if __name__ == "__main__":
    main(*sys.argv[1:2])
