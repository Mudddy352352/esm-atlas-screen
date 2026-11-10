#!/usr/bin/env python3
"""Fetch whichever NCBI searches have finished, and say what is still running.

Pairs with submit_all.py: every shortlisted protein already has a request id,
NCBI holds a finished result for 24 hours, and this collects them. Safe to run
as often as you like; a result already on disk is never fetched twice.

Usage: collect_nr.py <work_dir>
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_nr_blast import BLAST, TOOL  # noqa: E402


def call(params, timeout=300):
    req = urllib.request.Request(BLAST, data=urllib.parse.urlencode(params).encode())
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf8", "replace")


def status(rid):
    txt = call({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid,
                "TOOL": TOOL}, timeout=120)
    st = re.search(r"Status=(\w+)", txt)
    hits = "ThereAreHits=yes" in txt
    st = st.group(1) if st else "UNKNOWN"
    if st == "READY":
        return "ready" if hits else "no hits"
    return st.lower()


def result(rid):
    # Tabular is accepted and silently returns nothing; XML is the format that
    # actually carries the hits back.
    return call({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
                 "ALIGNMENTS": "100", "DESCRIPTIONS": "100", "TOOL": TOOL})


def main(work):
    res = os.path.join(work, "results")
    path = os.path.join(res, "state.json")
    state = json.load(open(path))
    waiting = []
    for acc, rec in sorted(state.items()):
        if rec.get("done"):
            continue
        rid = rec.get("rid")
        if not rid:
            continue
        st = status(rid)
        if st == "waiting":
            waiting.append(acc)
            print("[%s] still running" % acc, flush=True)
            continue
        if st == "ready":
            open(os.path.join(res, "%s.xml" % acc), "w").write(result(rid))
        rec.update(done=True, status=st,
                   collected=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        state[acc] = rec
        json.dump(state, open(path, "w"), indent=1)
        print("[%s] %s" % (acc, st), flush=True)
        time.sleep(3)
    done = sum(1 for v in state.values() if v.get("done"))
    print("collected %d of %d, still running: %d" % (done, len(state), len(waiting)))


if __name__ == "__main__":
    main(*sys.argv[1:2])
