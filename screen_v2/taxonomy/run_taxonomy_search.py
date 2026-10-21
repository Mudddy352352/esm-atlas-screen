#!/usr/bin/env python3
"""Step 6: search each tier-1 candidate structure against the taxonomy-bearing
Foldseek databases, so each candidate can be given a nearest labelled relative.

The MGnify protein database carries no taxonomic assignment for MGYP accessions
(the 2022_05 release ships no taxonomy file, and the Foldseek mgnify_esm30
database reports taxonomy=False), so the organism a candidate actually came from
is not recoverable. What this script produces instead is the taxonomy of the
closest structurally characterised relative, which is an inference about what the
protein resembles, not a claim about its source organism.

Queries are the real ESM Atlas structures, not ProstT5-predicted 3Di, so this is
a structure-vs-structure search.

Usage: run_taxonomy_search.py <outdir> <structure.pdb> [structure.pdb ...]
"""
import json
import os
import sys
import time
import urllib.error
import urllib.request

SERVER = "https://search.foldseek.com"
DATABASES = ["afdb50", "afdb-swissprot", "pdb100"]
MODE = "3diaa"


def post_multipart(url, fields, files):
    boundary = "----FoldseekTaxBoundary3b81d2"
    body = b""
    for k, v in fields:
        body += (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n").encode()
    for k, filename, content in files:
        body += (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"; "
                 f"filename=\"{filename}\"\r\nContent-Type: application/octet-stream\r\n\r\n").encode()
        body += content + b"\r\n"
    body += f"--{boundary}--\r\n".encode()
    req = urllib.request.Request(
        url, data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)


def get_json(url, timeout=300):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)


def submit(path):
    content = open(path, "rb").read()
    fields = [("mode", MODE)] + [("database[]", d) for d in DATABASES]
    return post_multipart(f"{SERVER}/api/ticket", fields,
                          [("q", os.path.basename(path), content)])


def wait(ticket_id, poll=15, limit=3600):
    waited = 0
    while waited < limit:
        try:
            st = get_json(f"{SERVER}/api/ticket/{ticket_id}", timeout=60).get("status")
        except Exception as exc:                     # transient server/proxy hiccup
            print(f"   poll error: {exc}", flush=True)
            st = None
        if st in ("COMPLETE", "ERROR"):
            return st
        time.sleep(poll)
        waited += poll
    return "TIMEOUT"


def main():
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    for path in sys.argv[2:]:
        name = os.path.basename(path).replace(".pdb", "")
        out = os.path.join(outdir, name + ".json")
        # completed results may have been gzipped after the run
        done = next((p for p in (out, out + ".gz")
                     if os.path.exists(p) and os.path.getsize(p) > 1000), None)
        if done:
            print(f"[skip] {name}", flush=True)
            continue
        try:
            ticket = submit(path)
        except Exception as exc:
            print(f"[fail] {name} submit: {exc}", flush=True)
            continue
        tid = ticket.get("id")
        print(f"[sent] {name} ticket={tid}", flush=True)
        status = wait(tid)
        if status != "COMPLETE":
            print(f"[fail] {name} status={status}", flush=True)
            continue
        try:
            res = get_json(f"{SERVER}/api/result/{tid}/0", timeout=900)
        except Exception as exc:
            print(f"[fail] {name} fetch: {exc}", flush=True)
            continue
        res["_ticket"] = tid
        res["_query_file"] = path
        res["_databases"] = DATABASES
        res["_mode"] = MODE
        res["_fetched_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with open(out, "w") as fh:
            json.dump(res, fh)
        counts = {b.get("db"): len(b.get("alignments") or []) for b in res.get("results", [])}
        print(f"[done] {name} {counts}", flush=True)


if __name__ == "__main__":
    main()
