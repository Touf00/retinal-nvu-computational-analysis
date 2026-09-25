import csv, gzip, json, glob
from pathlib import Path
from collections import defaultdict

art=Path("artifacts")
out=Path("combined"); out.mkdir(exist_ok=True)

summaries=[]
for fn in sorted(art.glob("chunk_*_summary.json")):
    summaries += json.loads(fn.read_text())

assert len(summaries)==31, f"Expected 31 libraries, got {len(summaries)}"
assert len({x["lib_id"] for x in summaries})==31
for x in summaries:
    assert int(x["observed_umi"])==int(x["expected_umi"]), x["lib_id"]

spec_counts=defaultdict(lambda: defaultdict(int))
meta={}
for s in summaries:
    specimen=s["bio_specimen"]
    if specimen in meta:
        assert meta[specimen]["pair"]==s["pair"]
        assert meta[specimen]["progression"]==s["progression"]
    else:
        meta[specimen]={"sample":specimen,"pair":str(s["pair"]),"progression":s["progression"]}
    fp=art/s["output"]
    assert fp.exists(), f"Missing {fp}"
    with gzip.open(fp,"rt") as f:
        r=csv.DictReader(f)
        for row in r:
            spec_counts[specimen][row["gene"]]+=int(row["count"])

assert len(meta)==28, f"Expected 28 specimens, got {len(meta)}"
pairs=defaultdict(set)
for m in meta.values(): pairs[m["pair"]].add(m["progression"])
assert len(pairs)==14 and all(v=={"Primary","Recurrent"} for v in pairs.values())

def pair_key(p):
    try: return int(p)
    except: return p
samples=sorted(meta, key=lambda s:(pair_key(meta[s]["pair"]), 0 if meta[s]["progression"]=="Primary" else 1, s))
genes=sorted(set().union(*(d.keys() for d in spec_counts.values())))

with open(out/"sample_metadata.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["sample","pair","progression"]); w.writeheader()
    for s in samples: w.writerow(meta[s])

with gzip.open(out/"pseudobulk_counts.csv.gz","wt",newline="") as f:
    w=csv.writer(f); w.writerow(["gene"]+samples)
    for g in genes:
        w.writerow([g]+[spec_counts[s].get(g,0) for s in samples])

# Per-specimen hard QC against frozen per-library UMI totals
expected_spec=defaultdict(int)
cells_spec=defaultdict(int)
for s in summaries:
    expected_spec[s["bio_specimen"]]+=int(s["expected_umi"])
    cells_spec[s["bio_specimen"]]+=int(s["n_ec99"])

with open(out/"specimen_qc.csv","w",newline="") as f:
    fields=["sample","pair","progression","n_ec99","expected_umi","matrix_umi","genes_detected","pass"]
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
    for s in samples:
        obs=sum(spec_counts[s].values())
        exp=expected_spec[s]
        row={"sample":s,"pair":meta[s]["pair"],"progression":meta[s]["progression"],
             "n_ec99":cells_spec[s],"expected_umi":exp,"matrix_umi":obs,
             "genes_detected":sum(v>0 for v in spec_counts[s].values()),
             "pass":obs==exp}
        assert row["pass"], row
        w.writerow(row)

with open(out/"library_qc.csv","w",newline="") as f:
    fields=["lib_id","bio_specimen","pair","progression","n_ec99","expected_umi","observed_umi","genes_detected","nul_lines"]
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
    for s in sorted(summaries,key=lambda x:x["lib_id"]): w.writerow({k:s[k] for k in fields})

audit={
    "libraries":len(summaries),"specimens":len(samples),"pairs":len(pairs),
    "genes_union":len(genes),"total_ec99_cells":sum(cells_spec.values()),
    "total_umi":sum(expected_spec.values()),"all_library_umi_checks":True,
    "all_specimen_umi_checks":True
}
(out/"combine_audit.json").write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2))
