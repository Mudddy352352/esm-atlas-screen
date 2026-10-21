#!/usr/bin/env python3
"""Submit query enzymes to the Foldseek web server and save raw results.

Queries are amino-acid FASTA; the server predicts 3Di with ProstT5 and runs a
structural search. Databases searched:
  mgnify_esm30 - ESM Metagenomic Atlas (MGnify proteins, 30% clustered)
  afdb50       - AlphaFold DB UniProt50 (what is already in sequence databases)
  pdb100       - PDB (identity check on the query itself)

Usage: run_foldseek.py <outdir> <query.fasta> [query.fasta ...]
"""
import json
import os
import sys
import time
import urllib.request

SERVER = "https://search.foldseek.com"
DATABASES = ["mgnify_esm30", "afdb50", "pdb100"]
MODE = "3diaa"


def post_multipart(url, fields, files):
    boundary = "----FoldseekBoundary7a3f9c"
    body = b""
    for k, v in fields:
        body += (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n").encode()
    for k, filename, content in files:
        body += (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"; "
                 f"filename=\"{filename}\"\r\nContent-Type: application/octet-stream\r\n\r\n").encode()
        body += content + b"\r\n"
    body += f"--{boundary}--\r\n".encode()
    req = urllib.request.Request(url, data=body,
                                headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)


def get_json(url, timeout=300):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)


def submit(fasta_path):
    content = open(fasta_path, "rb").read()
    fields = [("mode", MODE)] + [("database[]", d) for d in DATABASES]
    return post_multipart(f"{SERVER}/api/ticket", fields,
                          [("q", os.path.basename(fasta_path), content)])


def wait(ticket_id, poll=15, limit=3600):
    waited = 0
    while waited < limit:
        st = get_json(f"{SERVER}/api/ticket/{ticket_id}", timeout=60).get("status")
        if st in ("COMPLETE", "ERROR"):
            return st
        time.sleep(poll)
        waited += poll
    return "TIMEOUT"


def main():
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    for path in sys.argv[2:]:
        name = os.path.basename(path).replace(".fasta", "")
        out = os.path.join(outdir, name + ".json")
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            print(f"[skip] {name} already done", flush=True)
            continue
        t0 = time.time()
        try:
            ticket = submit(path)
        except Exception as e:
            print(f"[fail] {name} submit: {e}", flush=True)
            continue
        tid = ticket["id"]
        print(f"[sent] {name} ticket={tid}", flush=True)
        status = wait(tid)
        if status != "COMPLETE":
            print(f"[fail] {name} status={status}", flush=True)
            continue
        try:
            res = get_json(f"{SERVER}/api/result/{tid}/0", timeout=600)
        except Exception as e:
            print(f"[fail] {name} fetch: {e}", flush=True)
            continue
        res["_ticket"] = tid
        res["_query_file"] = path
        res["_databases"] = DATABASES
        res["_mode"] = MODE
        res["_fetched_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with open(out, "w") as fh:
            json.dump(res, fh)
        counts = {}
        for r in res.get("results", []):
            al = r["alignments"]
            al = al[0] if al and isinstance(al[0], list) else al
            counts[r["db"]] = len(al)
        print(f"[done] {name} {int(time.time()-t0)}s {counts}", flush=True)


if __name__ == "__main__":
    main()
