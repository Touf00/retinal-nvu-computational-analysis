import csv, gzip, json, os, re, subprocess, sys
from pathlib import Path
from collections import defaultdict

chunk=int(sys.argv[1])
NCHUNKS=int(sys.argv[2])
root=Path(__file__).resolve().parent
calls_path=root/"EC_cell_calls_FROZEN.csv.gz"
cluster_path=root/"lineage_cells_clustered.csv.gz"
out=Path("out"); out.mkdir(exist_ok=True)

PRIMARY_PAIRS={"1","12","13","16","20","22","25","27","29","32","35","40","41","9"}
BROAD_CLUSTERS={10,11,17,20,21,25,26}
DEFS=["EC95","EC99_ge2","EC99_broad"]

cluster={}
with gzip.open(cluster_path,"rt",newline="") as f:
    for r in csv.DictReader(f):
        if str(r["pair"]) in PRIMARY_PAIRS:
            cluster[r["cell_id"]]=int(float(r["cluster"]))

rows=[]
with gzip.open(calls_path,"rt",newline="") as f:
    for r in csv.DictReader(f):
        if str(r["pair"]) not in PRIMARY_PAIRS:
            continue
        r["pair"]=str(r["pair"])
        r["total_umi"]=int(float(r["total_umi"]))
        r["EC95"]=r["EC95"].lower()=="true"
        r["EC99"]=r["EC99"].lower()=="true"
        r["EC99_ge2"]=r["EC99_ge2"].lower()=="true"
        r["barcode"]=r["cell_id"].split("|")[-1]
        r["cluster"]=cluster.get(r["cell_id"])
        rows.append(r)

def qualifies(r,d):
    if d=="EC95": return r["EC95"]
    if d=="EC99_ge2": return r["EC99_ge2"]
    if d=="EC99_broad": return r["EC99"] and r["cluster"] in BROAD_CLUSTERS
    raise ValueError(d)

contracts={}
selected={}
for d in DEFS:
    z=[r for r in rows if qualifies(r,d)]
    spec=defaultdict(lambda:{"n":0,"umi":0,"pair":None,"progression":None})
    for r in z:
        s=spec[r["bio_specimen"]]
        s["n"]+=1; s["umi"]+=r["total_umi"]; s["pair"]=r["pair"]; s["progression"]=r["progression"]
    passed={sp for sp,v in spec.items() if v["n"]>=3 and v["umi"]>=5000}
    bypair=defaultdict(set)
    for sp in passed:
        v=spec[sp]; bypair[v["pair"]].add(v["progression"])
    pairs=sorted([p for p,v in bypair.items() if v=={"Primary","Recurrent"}],key=int)
    kept=[r for r in z if r["pair"] in set(pairs)]
    selected[d]=kept
    contracts[d]={
        "definition":d,
        "source":"frozen pre-expression cell calls",
        "minimum_cells":3,
        "minimum_umi":5000,
        "pair_ids":pairs,
        "n_pairs":len(pairs),
        "n_cells":len(kept),
        "n_specimens":len({r["bio_specimen"] for r in kept}),
        "broad_clusters":sorted(BROAD_CLUSTERS) if d=="EC99_broad" else None
    }

print("CONTRACTS",json.dumps(contracts,indent=2))
assert contracts["EC95"]["n_pairs"]==14
assert contracts["EC99_ge2"]["n_pairs"]==6
assert contracts["EC99_broad"]["n_pairs"]==13

# cells grouped by library and definition
bylib=defaultdict(lambda:defaultdict(list))
for d,z in selected.items():
    for r in z:
        bylib[r["lib_id"]][d].append(r)

libs=sorted(bylib)
todo=[lib for i,lib in enumerate(libs) if i % NCHUNKS == chunk]
print("CHUNK",chunk,"of",NCHUNKS,"libraries",len(todo),"of",len(libs))

def canon(s): return re.sub(r"[^ACGT]","",s)[-16:]
summaries=[]

for lib_id in todo:
    gsm,raw_spec,batch=lib_id.split(":")
    stem=f"{gsm}_{raw_spec}" + ("" if batch=="main" else f"_{batch}")
    bucket=gsm[:7]+"nnn"
    base=f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{bucket}/{gsm}/suppl"
    files={k:f"{stem}_{s}" for k,s in {
        "barcodes":"barcodes.tsv.gz","features":"features.tsv.gz","matrix":"matrix.mtx.gz"}.items()}
    for fn in files.values():
        subprocess.run(["curl","-fL","--retry","5","--retry-delay","2","-o",fn,f"{base}/{fn}"],check=True)

    defs_here=bylib[lib_id]
    memberships=defaultdict(list)
    expected={}
    metadata={}
    for d,cells in defs_here.items():
        expected[d]=sum(r["total_umi"] for r in cells)
        metadata[d]=cells[0]
        for r in cells: memberships[r["barcode"]].append(d)

    targets=set(memberships)
    with gzip.open(files["barcodes"],"rt",errors="replace") as f:
        bcs=[canon(x.strip()) for x in f if x.strip()]
    col_to_bc={i+1:b for i,b in enumerate(bcs) if b in targets}
    assert set(col_to_bc.values())==targets, f"{lib_id}: missing target barcodes"

    with gzip.open(files["features"],"rt",errors="replace") as f:
        genes=[x.replace("\x00","").strip().split("\t")[0] for x in f if x.strip()]

    counts={d:defaultdict(int) for d in defs_here}
    observed={d:0 for d in defs_here}
    nul_lines=0
    with gzip.open(files["matrix"],"rb") as f:
        for line in f:
            if not line.startswith(b"%"): break
        for line in f:
            if b"\x00" in line:
                nul_lines+=1; line=line.replace(b"\x00",b"")
            a=line.split()
            if len(a)<3: continue
            rr=int(a[0]); cc=int(a[1])
            bc=col_to_bc.get(cc)
            if bc is None: continue
            val=int(float(a[2])); gene=genes[rr-1]
            for d in memberships[bc]:
                observed[d]+=val
                counts[d][gene]+=val

    for d in defs_here:
        assert observed[d]==expected[d], f"{d} {lib_id}: UMI {observed[d]} != {expected[d]}"
        safe=lib_id.replace(":","__")
        op=out/f"{d}__{safe}.csv.gz"
        with gzip.open(op,"wt",newline="") as f:
            w=csv.writer(f); w.writerow(["gene","count"])
            for g,v in sorted(counts[d].items()):
                if v: w.writerow([g,v])
        cells=defs_here[d]; m=metadata[d]
        summaries.append({
            "definition":d,"lib_id":lib_id,"bio_specimen":m["bio_specimen"],
            "pair":m["pair"],"progression":m["progression"],
            "n_cells":len(cells),"expected_umi":expected[d],"observed_umi":observed[d],
            "genes_detected":sum(v>0 for v in counts[d].values()),"nul_lines":nul_lines,
            "output":op.name
        })
        print("PASS",d,lib_id,"cells",len(cells),"UMI",observed[d])

    for fn in files.values():
        try: Path(fn).unlink()
        except FileNotFoundError: pass

(out/f"sensitivity_chunk_{chunk}_summary.json").write_text(json.dumps(summaries,indent=2))
(out/f"sensitivity_chunk_{chunk}_contracts.json").write_text(json.dumps(contracts,indent=2))
print("DONE",chunk,"summaries",len(summaries))
