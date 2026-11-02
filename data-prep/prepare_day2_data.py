#!/usr/bin/env python3
"""
prepare_day2_data.py
Build the Day 2 data files (episodes/data/) from the nf-core/taxprofiler results
of the five Sun et al. (2021) simulated gut samples.

Why not use the Taxpasta tables directly? In this run Taxpasta was called without
a taxonomy (no --taxpasta_taxonomy_dir), so its tables hold taxonomy IDs only, with
all ranks mixed, and Taxpasta 0.7 does not support sylph. This script reads each
profiler's own per-sample report instead (which carries names), keeps species-level
rows, and writes one table per profiler with the same layout:

    taxonomy_id   name   Gut_sample1 ... Gut_sample5

Values are kept in each tool's native unit (MetaPhlAn and sylph: percent;
Bracken and mOTUs: read counts), so the lesson can show the conversion in R.

Inputs
  --results   taxprofiler results folder (contains metaphlan/, bracken/, ...)
  --truth     folder with Sun et al. 5_gut_taxonomic_abd.txt and 5_gut_sequence_abd.txt
  --cami      folder with the CAMI-format truth profiles (shenwei356/sun2021-cami-profiles,
              taxonomic_abd/Gut_sample*.profile), used for the NCBI taxonomy IDs
  --gtdb-map  species_map_gtdb_r220.tsv (name-based NCBI -> GTDB R220 map)
  --out       output folder (episodes/data)

Only the Python 3 standard library is used.
"""
import argparse, csv, glob, os, re, sys
from collections import OrderedDict

SAMPLES = [f"Gut_sample{i}" for i in range(1, 6)]


def sample_of(path):
    m = re.search(r"(Gut_sample\d+)", os.path.basename(path))
    if not m:
        sys.exit("Cannot find a Gut_sampleN name in " + path)
    return m.group(1)


def one(pattern):
    hits = sorted(glob.glob(pattern))
    if len(hits) != len(SAMPLES):
        sys.exit("Expected %d files for %s, found %d" % (len(SAMPLES), pattern, len(hits)))
    return hits


def write_table(path, rows, value_fmt, extra_cols=()):
    """rows: OrderedDict id -> dict(name=..., <extra>=..., Gut_sampleN=value)"""
    cols = ["taxonomy_id", "name", *extra_cols, *SAMPLES]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for tid, r in rows.items():
            w.writerow([tid, r["name"], *[r.get(c, "") for c in extra_cols],
                        *[value_fmt(r.get(s, 0)) for s in SAMPLES]])
    print("  wrote %-28s %4d rows" % (os.path.basename(path), len(rows)), file=sys.stderr)


def by_total(rows):
    return OrderedDict(sorted(rows.items(), key=lambda kv: -sum(kv[1].get(s, 0) for s in SAMPLES)))


def pct(x):
    return "%.5f" % x if x else "0"


def count(x):
    return "%d" % round(x)


def frac(x):
    return "%.8f" % x if x else "0"


# ----------------------------------------------------------------------------- profilers
def metaphlan(results):
    rows = {}
    for f in one(f"{results}/metaphlan/*/*metaphlan_profile.txt"):
        s = sample_of(f)
        for line in open(f):
            if line.startswith("#"):
                continue
            clade, taxpath, value = line.rstrip("\n").split("\t")[:3]
            last = clade.split("|")[-1]
            if not last.startswith("s__"):
                continue
            name = last[3:].replace("_", " ")
            tid = taxpath.split("|")[-1] or "MetaPhlAn:" + last[3:]  # SGBs without an NCBI species
            rows.setdefault(tid, {"name": name})[s] = float(value)
    return by_total(rows)


def bracken(results):
    rows = {}
    for f in one(f"{results}/bracken/*/*.bracken.tsv"):
        s = sample_of(f)
        for r in csv.DictReader(open(f), delimiter="\t"):
            rows.setdefault(r["taxonomy_id"], {"name": r["name"]})[s] = int(r["new_est_reads"])
    return by_total(rows)


def motus(results):
    rows = {}
    for f in one(f"{results}/motus/*/*.out"):
        s = sample_of(f)
        for line in open(f):
            if line.startswith("#"):
                continue
            label, tid, value = line.rstrip("\n").split("\t")[:3]
            v = float(value)
            if v == 0:
                continue
            if label.startswith("unassigned"):
                tid, name = "-1", "unassigned"
            else:
                name = re.sub(r"\s*\[[^\]]*\]$", "", label)
                if tid in ("NA", ""):  # mOTUs without an NCBI species: keep the mOTU ID
                    tid = "mOTUs:" + re.search(r"\[([^\]]*)\]$", label).group(1)
            e = rows.setdefault(tid, {"name": name})
            e[s] = e.get(s, 0) + v
    return by_total(rows)


def sylph(results):
    tax, seq = {}, {}
    for f in one(f"{results}/sylph/*/*.sylphmpa"):
        s = sample_of(f)
        for line in open(f):
            if line.startswith(("#", "clade_name")):
                continue
            p = line.rstrip("\n").split("\t")
            last = p[0].split("|")[-1]
            if not last.startswith("s__"):
                continue
            sp = last[3:]
            for d, v in ((tax, p[1]), (seq, p[2])):
                d.setdefault(sp, {"name": sp, "gtdb_lineage": p[0].replace("|", ";")})[s] = float(v)
    return by_total(tax), by_total(seq)


def kraken2_unclassified(results):
    out = {}
    for f in one(f"{results}/kraken2/*/*kraken2.report.txt"):
        u = c = 0
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if p[4].strip() == "0":
                u = int(p[1])
            elif p[4].strip() == "1":
                c = int(p[1])
        out[sample_of(f)] = (u, u + c)
    return out


def read_pairs(results):
    import json
    out = {}
    for f in one(f"{results}/fastp/*.fastp.json"):
        d = json.load(open(f))["summary"]
        out[sample_of(f)] = (d["before_filtering"]["total_reads"] // 2, d["after_filtering"]["total_reads"] // 2)
    return out


# ----------------------------------------------------------------------------- truth
def read_sun(path):
    rows = list(csv.reader(open(path), delimiter="\t"))
    return OrderedDict((r[0], [float(x) for x in r[1:]]) for r in rows[1:] if r)


def cami_species(cami_dir):
    d = {}
    for f in one(f"{cami_dir}/Gut_sample*.profile"):
        s = sample_of(f)
        for line in open(f):
            if line.startswith("@") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[1] == "species":
                d.setdefault(p[0], {"name": p[3].split("|")[-1]})[s] = float(p[4]) / 100
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--cami", required=True)
    ap.add_argument("--gtdb-map", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    print("Profilers:", file=sys.stderr)
    mpa, brk, mot = metaphlan(a.results), bracken(a.results), motus(a.results)
    syl_tax, syl_seq = sylph(a.results)
    write_table(f"{a.out}/metaphlan.tsv", mpa, pct)
    write_table(f"{a.out}/motus.tsv", mot, count)
    write_table(f"{a.out}/bracken.tsv", brk, count)
    write_table(f"{a.out}/sylph_taxonomic.tsv", syl_tax, pct, extra_cols=("gtdb_lineage",))
    write_table(f"{a.out}/sylph_sequence.tsv", syl_seq, pct, extra_cols=("gtdb_lineage",))

    # truth: Sun names -> NCBI taxonomy IDs (CAMI profiles), matched on the abundance vectors
    tax = read_sun(f"{a.truth}/5_gut_taxonomic_abd.txt")
    seq = read_sun(f"{a.truth}/5_gut_sequence_abd.txt")
    cami = cami_species(a.cami)
    sun2tid = {}
    for sp, v in tax.items():
        tid = min(cami, key=lambda t: sum(abs(v[i] - cami[t].get(SAMPLES[i], 0)) for i in range(5)))
        err = sum(abs(v[i] - cami[tid].get(SAMPLES[i], 0)) for i in range(5))
        if err > 1e-6:
            sys.exit("No exact CAMI match for %s (closest %s, L1 %.2g)" % (sp, cami[tid]["name"], err))
        sun2tid[sp] = tid
    if len(set(sun2tid.values())) != len(sun2tid):
        sys.exit("Two truth species matched the same taxonomy ID")

    # GTDB R220 names: name-based map, then genome-level evidence from sylph for split species.
    gmap = {r["truth_species"]: r["gtdb_lineage"].split(";")[-1][3:]
            for r in csv.DictReader(open(a.gtdb_map), delimiter="\t")}
    gmap_by_name = dict(gmap)
    evidence = {sp: "name (NCBI type strain in GTDB R220)" for sp in tax}
    unmatched_truth = [sp for sp in tax if gmap[sp] not in syl_tax]
    unclaimed = [g for g in syl_tax if g not in set(gmap.values())]
    for sp in unmatched_truth:
        v = tax[sp]
        def dist(g):
            return sum(abs(v[i] * 100 - syl_tax[g].get(SAMPLES[i], 0)) for i in range(5))
        g = min(unclaimed, key=dist)
        same_presence = all((v[i] > 0) == (syl_tax[g].get(SAMPLES[i], 0) > 0) for i in range(5))
        if not (same_presence and dist(g) < 0.5):  # percent points summed over 5 samples
            print("  WARNING: no clear GTDB match for %s; keeping %s" % (sp, gmap[sp]), file=sys.stderr)
            continue
        print("  %-30s GTDB R220 %-32s -> genome in %s" % (sp, gmap[sp], g), file=sys.stderr)
        gmap[sp] = g
        unclaimed.remove(g)
        evidence[sp] = "genome (sylph ANI and abundance pattern)"
    for sp in tax:
        if evidence[sp].startswith("name") and gmap[sp] in syl_tax:
            evidence[sp] = "name, confirmed by sylph detection"

    ncbi_name = {sp: cami[sun2tid[sp]]["name"] for sp in tax}
    for kind, data in (("taxonomic", tax), ("sequence", seq)):
        rows = OrderedDict()
        for sp, v in sorted(data.items(), key=lambda kv: -sum(kv[1])):
            rows[sun2tid[sp]] = {"name": ncbi_name[sp], "gtdb_r220_species": gmap[sp],
                                 **{SAMPLES[i]: v[i] for i in range(5)}}
        write_table(f"{a.out}/truth_{kind}.tsv", rows, frac, extra_cols=("gtdb_r220_species",))

    with open(f"{a.out}/ncbi_to_gtdb_r220.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["taxonomy_id", "ncbi_name", "sun2021_name", "gtdb_r220_species",
                    "gtdb_r220_species_by_name", "evidence"])
        for sp in sorted(tax, key=lambda s: ncbi_name[s]):
            w.writerow([sun2tid[sp], ncbi_name[sp], sp, gmap[sp], gmap_by_name[sp], evidence[sp]])
    print("  wrote ncbi_to_gtdb_r220.tsv       %4d rows" % len(tax), file=sys.stderr)

    k2, pairs = kraken2_unclassified(a.results), read_pairs(a.results)
    with open(f"{a.out}/metadata.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "community", "n_species_truth", "read_pairs_raw", "read_pairs_qc",
                    "kraken2_unclassified_pct"])
        for i, s in enumerate(SAMPLES):
            n = sum(1 for v in tax.values() if v[i] > 0)
            w.writerow([s, "simulated human gut", n, pairs[s][0], pairs[s][1],
                        "%.2f" % (100 * k2[s][0] / k2[s][1])])
    print("  wrote metadata.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
