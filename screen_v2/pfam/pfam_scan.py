#!/usr/bin/env python3
"""Pfam-profile comparator for the Foldseek screen.

Question: would a Pfam profile search have connected each structure-only
candidate to the same protein family as the drug-metabolising query enzyme?

Run hmmscan (pyhmmer) at Pfam gathering thresholds over
  - the 9 full-length query enzymes
  - the 1,640 distinct structure-only candidates
then ask, per candidate, whether it shares at least one Pfam family with its
own query.

Outputs (written to the directory given as --out):
  pfam_query_domains.tsv      one row per query x Pfam domain
  pfam_candidate_domains.tsv  one row per candidate x Pfam domain
  pfam_recovery.tsv           one row per candidate x query pair, with verdict
  pfam_summary.txt            the counts to quote in the manuscript
"""
import argparse, collections, csv, gzip, os, sys

import pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet
from pyhmmer.plan7 import HMMFile


def asstr(x):
    """pyhmmer returns names as bytes in some versions and str in others."""
    if x is None:
        return ""
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def load_seqs(path, alphabet):
    with SequenceFile(path, digital=True, alphabet=alphabet) as fh:
        return list(fh)


def scan(seqs, hmm_path, cpus):
    """Return {seq_name: {pfam_name: (accession, evalue, bitscore)}}."""
    out = collections.defaultdict(dict)
    with HMMFile(hmm_path) as hf:
        hmms = list(hf)
    print("loaded %d Pfam models" % len(hmms), flush=True)
    for hits in pyhmmer.hmmer.hmmscan(seqs, hmms, cpus=cpus,
                                      bit_cutoffs="gathering"):
        qname = asstr(hits.query.name)
        for hit in hits:
            if not hit.included:
                continue
            name = asstr(hit.name)
            acc = asstr(hit.accession)
            prev = out[qname].get(name)
            if prev is None or hit.evalue < prev[1]:
                out[qname][name] = (acc, hit.evalue, hit.score)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--pairs", required=True,
                    help="pass-1 candidates.tsv (query/target pairs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpus", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    alpha = Alphabet.amino()
    queries = load_seqs(args.queries, alpha)
    cands = load_seqs(args.candidates, alpha)
    print("queries %d, candidates %d" % (len(queries), len(cands)), flush=True)

    qdom = scan(queries, args.pfam, args.cpus)
    print("scanned queries", flush=True)
    cdom = scan(cands, args.pfam, args.cpus)
    print("scanned candidates", flush=True)

    def write_domains(path, table):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sequence", "pfam_name", "pfam_acc", "evalue", "bits"])
            for name in sorted(table):
                for dom, (acc, ev, bits) in sorted(table[name].items(),
                                                   key=lambda kv: kv[1][1]):
                    w.writerow([name, dom, acc, "%.3g" % ev, "%.1f" % bits])

    write_domains(os.path.join(args.out, "pfam_query_domains.tsv"), qdom)
    write_domains(os.path.join(args.out, "pfam_candidate_domains.tsv"), cdom)

    # query model name in the screen (e.g. P05804_uidA_gus_w214-603) maps to the
    # full-length query sequence name by its accession prefix.
    acc_of_model = {}
    pairs = []
    with open(args.pairs) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            model = row["query"]
            acc = model.split("_")[0]
            acc_of_model[model] = acc
            pairs.append((model, row["target"].replace(".pdb.gz", "")))

    qdom_by_acc = {}
    for qname, doms in qdom.items():
        qdom_by_acc[qname.split("_")[0]] = set(doms)

    rows = []
    for model, target in pairs:
        acc = acc_of_model[model]
        qset = qdom_by_acc.get(acc, set())
        cset = set(cdom.get(target, {}))
        shared = sorted(qset & cset)
        rows.append({
            "query_model": model,
            "query_acc": acc,
            "candidate": target,
            "query_pfam": ";".join(sorted(qset)) or "-",
            "candidate_pfam": ";".join(sorted(cset)) or "-",
            "shared_pfam": ";".join(shared) or "-",
            "candidate_has_pfam": "yes" if cset else "no",
            "pfam_recovers_link": "yes" if shared else "no",
        })

    with open(os.path.join(args.out, "pfam_recovery.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    n_pairs = len(rows)
    n_recov = sum(r["pfam_recovers_link"] == "yes" for r in rows)
    n_haspfam = sum(r["candidate_has_pfam"] == "yes" for r in rows)
    prot = {}
    for r in rows:
        p = prot.setdefault(r["candidate"], {"rec": False, "pfam": False})
        p["rec"] |= r["pfam_recovers_link"] == "yes"
        p["pfam"] |= r["candidate_has_pfam"] == "yes"

    per_query = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        per_query[r["query_model"]][0] += 1
        per_query[r["query_model"]][1] += r["pfam_recovers_link"] == "yes"

    noq = sorted({r["query_acc"] for r in rows if r["query_pfam"] == "-"})

    lines = []
    lines.append("Pfam profile comparator, hmmscan at gathering thresholds")
    lines.append("Pfam file: %s" % os.path.basename(args.pfam))
    lines.append("")
    lines.append("candidate (query, protein) pairs: %d" % n_pairs)
    lines.append("distinct candidate proteins: %d" % len(prot))
    lines.append("pairs where the candidate has any Pfam domain: %d (%.1f%%)"
                 % (n_haspfam, 100.0 * n_haspfam / n_pairs))
    lines.append("pairs where candidate shares a Pfam family with its query: "
                 "%d (%.1f%%)" % (n_recov, 100.0 * n_recov / n_pairs))
    lines.append("distinct proteins sharing a Pfam family with any query: "
                 "%d (%.1f%%)" % (sum(p["rec"] for p in prot.values()),
                                  100.0 * sum(p["rec"] for p in prot.values()) / len(prot)))
    lines.append("distinct proteins with any Pfam domain: %d (%.1f%%)"
                 % (sum(p["pfam"] for p in prot.values()),
                    100.0 * sum(p["pfam"] for p in prot.values()) / len(prot)))
    if noq:
        lines.append("queries with NO Pfam domain at gathering threshold: %s"
                     % ", ".join(noq))
    lines.append("")
    lines.append("per query model: candidates, recovered by Pfam, percent")
    for m in sorted(per_query):
        tot, rec = per_query[m]
        lines.append("  %-32s %5d %5d %6.1f%%" % (m, tot, rec, 100.0 * rec / tot))
    lines.append("")
    lines.append("query Pfam domains:")
    for qname in sorted(qdom):
        lines.append("  %-32s %s" % (qname, ";".join(sorted(qdom[qname])) or "-"))

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(args.out, "pfam_summary.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
