import csv, gzip, json, os, re, subprocess, sys
from pathlib import Path
from collections import defaultdict

chunk=int(sys.argv[1])
root=Path(__file__).resolve().parent
cfg=json.loads((root/f"input_chunk_{chunk}.json").read_text())
out=Path("out"); out.mkdir(exist_ok=True)

def canon(s):
    return re.sub(r"[^ACGT]","",s)[-16:]

summaries=[]
for lib in cfg["libraries"]:
    gsm=lib["gsm"]; spec=lib["raw_spec"]; batch=lib["batch"]
    stem=f"{gsm}_{spec}" + ("" if batch=="main" else f"_{batch}")
    bucket=gsm[:7]+"nnn"
    base=f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{bucket}/{gsm}/suppl"
    files={
        "barcodes":f"{stem}_barcodes.tsv.gz",
        "features":f"{stem}_features.tsv.gz",
        "matrix":f"{stem}_matrix.mtx.gz",
    }
    for key,fn in files.items():
        if not Path(fn).exists():
            subprocess.run(["curl","-fL","--retry","5","--retry-delay","2","-o",fn,f"{base}/{fn}"],check=True)

    targets=set(lib["barcodes"])
    with gzip.open(files["barcodes"],"rt",errors="replace") as f:
        raw=[x.strip() for x in f if x.strip()]
    bcs=[canon(x) for x in raw]
    col_to_target={i+1:b for i,b in enumerate(bcs) if b in targets}
    found=set(col_to_target.values())
    assert found==targets, f"{lib['lib_id']}: missing frozen barcodes {sorted(targets-found)}"

    with gzip.open(files["features"],"rt",errors="replace") as f:
        genes=[x.replace("\x00","").strip().split("\t")[0] for x in f if x.strip()]

    counts=defaultdict(int); total=0; nul_lines=0
    with gzip.open(files["matrix"],"rb") as f:
        for line in f:
            if not line.startswith(b"%"):
                break
        for line in f:
            if b"\x00" in line:
                nul_lines+=1
                line=line.replace(b"\x00",b"")
            a=line.split()
            if len(a)<3: continue
            rr=int(a[0]); cc=int(a[1])
            if cc not in col_to_target: continue
            val=int(float(a[2]))
            total+=val
            counts[genes[rr-1]]+=val

    expected=int(lib["expected_umi"])
    assert total==expected, f"{lib['lib_id']}: UMI mismatch observed={total} expected={expected}"
    assert len(targets)==int(lib["n_ec99"])

    safe=lib["lib_id"].replace(":","__")
    op=out/f"{safe}.csv.gz"
    with gzip.open(op,"wt",newline="") as f:
        w=csv.writer(f); w.writerow(["gene","count"])
        for g,v in sorted(counts.items()):
            if v: w.writerow([g,v])

    s={k:lib[k] for k in ["lib_id","bio_specimen","pair","progression","n_ec99","expected_umi"]}
    s.update({"observed_umi":total,"genes_detected":sum(v>0 for v in counts.values()),
              "nul_lines":nul_lines,"output":op.name})
    summaries.append(s)
    print("PASS",lib["lib_id"],"cells",len(targets),"UMI",total,"genes",s["genes_detected"],"NUL",nul_lines)

    for fn in files.values():
        try: Path(fn).unlink()
        except FileNotFoundError: pass

(out/f"chunk_{chunk}_summary.json").write_text(json.dumps(summaries,indent=2))
print("CHUNK COMPLETE",chunk,"libraries",len(summaries))
