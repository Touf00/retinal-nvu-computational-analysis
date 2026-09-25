import csv, gzip, json
from pathlib import Path
from collections import defaultdict

art=Path("artifacts")
root=Path("sensitivity"); root.mkdir(exist_ok=True)

summaries=[]
contract_sets=[]
for fn in sorted(art.glob("sensitivity_chunk_*_summary.json")):
    summaries += json.loads(fn.read_text())
for fn in sorted(art.glob("sensitivity_chunk_*_contracts.json")):
    contract_sets.append(json.loads(fn.read_text()))
assert contract_sets, "No contracts"
contracts=contract_sets[0]
assert all(x==contracts for x in contract_sets), "Contract mismatch across chunks"

defs=sorted(contracts)
for d in defs:
    out=root/d; out.mkdir(parents=True,exist_ok=True)
    ss=[x for x in summaries if x["definition"]==d]
    assert all(int(x["expected_umi"])==int(x["observed_umi"]) for x in ss)
    spec_counts=defaultdict(lambda:defaultdict(int))
    meta={}
    for s in ss:
        sp=s["bio_specimen"]
        meta.setdefault(sp,{"sample":sp,"pair":str(s["pair"]),"progression":s["progression"]})
        fp=art/s["output"]
        assert fp.exists(), fp
        with gzip.open(fp,"rt") as f:
            for r in csv.DictReader(f):
                spec_counts[sp][r["gene"]]+=int(r["count"])
    pairs=contracts[d]["pair_ids"]
    keep_samples=[sp for sp,m in meta.items() if str(m["pair"]) in set(map(str,pairs))]
    assert len(keep_samples)==2*len(pairs)
    def pk(s):
        m=meta[s]
        return (int(m["pair"]),0 if m["progression"]=="Primary" else 1,s)
    samples=sorted(keep_samples,key=pk)
    genes=sorted(set().union(*(spec_counts[s].keys() for s in samples)))
    with open(out/"sample_metadata.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["sample","pair","progression"]); w.writeheader()
        for s in samples: w.writerow(meta[s])
    with gzip.open(out/"pseudobulk_counts.csv.gz","wt",newline="") as f:
        w=csv.writer(f); w.writerow(["gene"]+samples)
        for g in genes: w.writerow([g]+[spec_counts[s].get(g,0) for s in samples])
    expected_spec=defaultdict(int); cells_spec=defaultdict(int)
    for s in ss:
        if str(s["pair"]) in set(map(str,pairs)):
            expected_spec[s["bio_specimen"]]+=int(s["expected_umi"])
            cells_spec[s["bio_specimen"]]+=int(s["n_cells"])
    with open(out/"specimen_qc.csv","w",newline="") as f:
        fields=["sample","pair","progression","n_cells","expected_umi","matrix_umi","genes_detected","pass"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for s in samples:
            obs=sum(spec_counts[s].values()); exp=expected_spec[s]
            row={"sample":s,"pair":meta[s]["pair"],"progression":meta[s]["progression"],
                 "n_cells":cells_spec[s],"expected_umi":exp,"matrix_umi":obs,
                 "genes_detected":sum(v>0 for v in spec_counts[s].values()),"pass":obs==exp}
            assert row["pass"],row; w.writerow(row)
    (out/"contract.json").write_text(json.dumps(contracts[d],indent=2))
    print(d,"pairs",len(pairs),"samples",len(samples),"genes",len(genes),
          "cells",sum(cells_spec.values()),"UMI",sum(expected_spec.values()))

(root/"contracts.json").write_text(json.dumps(contracts,indent=2))
