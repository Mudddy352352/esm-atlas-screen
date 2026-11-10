#!/usr/bin/env python3
"""Does a shared fold carry the active site? The same question the screen asks,
asked across many enzyme families instead of nine.

For each sampled M-CSA entry (curated catalytic residues, reference PDB) this
searches AlphaFold DB with Foldseek, then measures for every same-fold hit
whether the curated catalytic residues are still there and still correctly
placed, splitting hits into the ones a single-sequence search also finds and
the ones only structure search reaches. That is the paper's 83.9% against
42.3% comparison, generalised.

Resumable: one row per hit appended to the results TSV, one line per finished
entry in done.txt, and the Foldseek JSON is deleted as soon as it is scored so
the scan stays inside the disk allowance.

Usage: mcsa_retention_scan.py <work_dir> <n_entries> [max_minutes]
"""
import json
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd
import pyhmmer

sys.path.insert(0, _ROOT + "/screen_v2/code")
from active_site_check import ca_from_pdb, best_fit, walk  # noqa: E402

# Repository root. Every path below is relative to it, so the release works
# wherever it is cloned. Override with ESMATLAS_ROOT to point at another copy.
_ROOT = os.environ.get("ESMATLAS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

MCSA = "https://www.ebi.ac.uk/thornton-srv/m-csa/api"
FOLDSEEK = "https://search.foldseek.com/api"
RCSB = "https://files.rcsb.org/download/%s.pdb"
DB = "afdb50"
RMSD_OK, MAD_OK, TM_MIN, PHMMER_E = 2.0, 1.5, 0.5, 1e-3
THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
         "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
         "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
         "Tyr": "Y", "Val": "V"}


def fetch(url, tries=3, timeout=120):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh:
                return fh.read()
        except Exception as exc:
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))
            del exc


def mcsa_entries(limit):
    out, page = [], 1
    while len(out) < limit * 4 and page <= 12:
        d = json.loads(fetch("%s/entries/?format=json&page=%d&page_size=100" % (MCSA, page)))
        out.extend(d["results"])
        if not d.get("next"):
            break
        page += 1
    return out


def catalytic_residues(entry):
    """Reference PDB chain and its curated catalytic residues."""
    best = {}
    for r in entry.get("residues", []):
        for rc in r.get("residue_chains", []):
            if not rc.get("is_reference"):
                continue
            key = (rc["pdb_id"], rc["chain_name"])
            best.setdefault(key, []).append((int(rc["auth_resid"]),
                                             THREE.get(rc["code"], "X")))
    if not best:
        return None, None, []
    (pdb, chain), res = max(best.items(), key=lambda kv: len(kv[1]))
    res = sorted(set(res))
    return pdb, chain, [(p, a) for p, a in res if a != "X"]


def write_chain(pdb_id, chain, path):
    raw = fetch(RCSB % pdb_id.upper()).decode("utf-8", "replace")
    keep, seen = [], set()
    for line in raw.splitlines():
        if line.startswith("ATOM") and line[21] == chain:
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            keep.append(line)
            seen.add(line[22:27].strip())
    if len(keep) < 100:
        return False
    with open(path, "w") as fh:
        fh.write("\n".join(keep) + "\nEND\n")
    return True


def ca_index(path, chain):
    """auth residue number -> (0-based CA index, one-letter code)."""
    idx, i = {}, 0
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            num = line[22:27].strip()
            aa = THREE.get(line[17:20].strip().capitalize(), "X")
            idx[num] = (i, aa)
            i += 1
    return idx


def foldseek(path):
    import mimetypes  # noqa: F401
    boundary = "----mcsa%d" % int(time.time() * 1000)
    parts = []
    with open(path, "rb") as fh:
        blob = fh.read()
    parts.append(("q", os.path.basename(path), blob))
    body = b""
    for name, fname, data in parts:
        body += ("--%s\r\nContent-Disposition: form-data; name=\"%s\"; filename=\"%s\"\r\n"
                 "Content-Type: application/octet-stream\r\n\r\n" % (boundary, name, fname)).encode()
        body += data + b"\r\n"
    for key, val in [("mode", "3diaa"), ("database[]", DB)]:
        body += ("--%s\r\nContent-Disposition: form-data; name=\"%s\"\r\n\r\n%s\r\n"
                 % (boundary, key, val)).encode()
    body += ("--%s--\r\n" % boundary).encode()
    req = urllib.request.Request(FOLDSEEK + "/ticket", data=body, method="POST")
    req.add_header("Content-Type", "multipart/form-data; boundary=" + boundary)
    with urllib.request.urlopen(req, timeout=180) as fh:
        ticket = json.load(fh)
    tid = ticket["id"]
    for _ in range(240):
        st = json.loads(fetch("%s/ticket/%s" % (FOLDSEEK, tid), timeout=60))
        if st.get("status") == "COMPLETE":
            break
        if st.get("status") == "ERROR":
            raise RuntimeError("foldseek error for " + path)
        time.sleep(5)
    else:
        raise TimeoutError("foldseek timeout for " + path)
    return json.loads(fetch("%s/result/%s/0" % (FOLDSEEK, tid), timeout=600))


def phmmer_blind(qseq, qname, targets):
    alpha = pyhmmer.easel.Alphabet.amino()
    q = pyhmmer.easel.TextSequence(name=qname.encode(), sequence=qseq).digitize(alpha)
    seqs = [pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(alpha)
            for k, v in targets.items() if v]
    if not seqs:
        return {}
    db = pyhmmer.easel.DigitalSequenceBlock(alpha, seqs)
    found = {}
    for top in pyhmmer.hmmer.phmmer([q], db, E=10.0, cpus=0):
        for hit in top:
            name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
            found[name] = hit.evalue
    return found


def score_entry(entry, work, rows_path):
    mcsa_id = entry["mcsa_id"]
    pdb, chain, res = catalytic_residues(entry)
    if not pdb or len(res) < 2:
        return 0, "no reference chain with two or more catalytic residues"
    qpath = os.path.join(work, "queries", "mcsa%04d_%s_%s.pdb" % (mcsa_id, pdb, chain))
    if not os.path.exists(qpath) and not write_chain(pdb, chain, qpath):
        return 0, "chain too short or not downloadable"
    idx = ca_index(qpath, chain)
    wpos, codes = [], []
    for num, aa in res:
        hit = idx.get(str(num))
        if hit is None or hit[1] != aa:
            continue
        wpos.append(hit[0] + 1)
        codes.append(aa)
    if len(wpos) < 2:
        return 0, "catalytic residues not found in the chain"

    qxyz, qseq = ca_from_pdb(qpath)
    order = np.argsort(wpos)
    wpos = list(np.array(wpos)[order])
    codes = list(np.array(codes)[order])
    qs = qxyz[np.array(wpos) - 1]
    qd = np.linalg.norm(qs[:, None] - qs[None, :], axis=-1)

    res_json = foldseek(qpath)
    al = []
    for block in res_json.get("results", []):
        if block.get("db") != DB:
            continue
        a = block["alignments"]
        al = a[0] if a and isinstance(a[0], list) else a
    if not al:
        return 0, "no hits"

    targets = {str(a["target"]).split()[0]: (a.get("tSeq") or "") for a in al}
    ph = phmmer_blind(qseq, "mcsa%d" % mcsa_id, targets)

    out = []
    for a in al:
        try:
            txyz = np.array([float(v) for v in a["tCa"].split(",")]).reshape(-1, 3)
        except Exception:
            continue
        pairs = walk(a)
        qp = np.array([qxyz[i - 1] for i in pairs if i <= len(qxyz)])
        tp = np.array([txyz[pairs[i][0] - 1] for i in pairs
                       if i <= len(qxyz) and pairs[i][0] <= len(txyz)])
        if len(qp) < 3 or len(qp) != len(tp):
            continue
        R, t, tm = best_fit(qp, tp, a["qLen"])
        if R is None or tm < TM_MIN:
            continue
        n_id = n_al = 0
        devs, tsite, keep = [], [], []
        for k, (p, aa) in enumerate(zip(wpos, codes)):
            hit = pairs.get(int(p))
            if hit is None or hit[0] > len(txyz):
                continue
            tpos, tres = hit
            n_al += 1
            n_id += (tres == aa)
            devs.append(float(np.linalg.norm(R @ txyz[tpos - 1] + t - qxyz[int(p) - 1])))
            tsite.append(txyz[tpos - 1])
            keep.append(k)
        mad = np.nan
        if len(keep) >= 3:
            ts = np.array(tsite)
            td = np.linalg.norm(ts[:, None] - ts[None, :], axis=-1)
            sub = np.ix_(keep, keep)
            mad = float(np.mean(np.abs(qd[sub][np.triu_indices(len(keep), 1)]
                                       - td[np.triu_indices(len(keep), 1)])))
        rmsd = float(np.sqrt(np.mean(np.square(devs)))) if devs else np.nan
        tid = str(a["target"]).split()[0]
        ev = ph.get(tid)
        out.append(dict(
            mcsa_id=mcsa_id, enzyme=entry.get("enzyme_name", ""),
            pdb=pdb, chain=chain, n_catalytic=len(wpos),
            target=tid, target_desc=" ".join(str(a["target"]).split()[1:]),
            seq_id_pct=round(float(a.get("seqId", 0)), 1),
            tm_score_qnorm=round(float(tm), 4),
            n_cat_aligned=n_al, n_cat_identical=n_id,
            site_rmsd=round(rmsd, 2) if rmsd == rmsd else np.nan,
            site_dist_mad=round(mad, 2) if mad == mad else np.nan,
            phmmer_evalue=ev,
            seq_detectable=bool(ev is not None and ev <= PHMMER_E),
            catalytic_kept=bool(n_al == len(wpos) and n_id == len(wpos)
                                and rmsd == rmsd and rmsd <= RMSD_OK
                                and (mad != mad or mad <= MAD_OK))))
    if not out:
        return 0, "no scorable hits"
    df = pd.DataFrame(out)
    df.to_csv(rows_path, sep="\t", index=False, mode="a",
              header=not os.path.exists(rows_path))
    return len(df), "ok"


def main(work, n_entries, max_minutes=240):
    n_entries, max_minutes = int(n_entries), float(max_minutes)
    os.makedirs(os.path.join(work, "queries"), exist_ok=True)
    os.makedirs(os.path.join(work, "results"), exist_ok=True)
    rows_path = os.path.join(work, "results", "mcsa_hits.tsv")
    done_path = os.path.join(work, "results", "done.txt")
    done = set()
    if os.path.exists(done_path):
        done = {int(x.split("\t")[0]) for x in open(done_path) if x.strip()}

    entries = mcsa_entries(n_entries)
    print("M-CSA entries listed:", len(entries), flush=True)
    started = time.time()
    n_ok = 0
    for e in entries:
        if n_ok >= n_entries or (time.time() - started) / 60 > max_minutes:
            break
        if e["mcsa_id"] in done:
            continue
        try:
            full = json.loads(fetch("%s/entries/%d/?format=json" % (MCSA, e["mcsa_id"])))
            n, why = score_entry(full, work, rows_path)
        except Exception as exc:
            n, why = 0, "failed: %s" % exc
        with open(done_path, "a") as fh:
            fh.write("%d\t%d\t%s\n" % (e["mcsa_id"], n, why))
        print("[mcsa %d] %s: %d hits (%s)" % (e["mcsa_id"], e.get("enzyme_name", "")[:40], n, why),
              flush=True)
        if n:
            n_ok += 1
        time.sleep(3)
    print("finished, entries scored this run:", n_ok, flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
