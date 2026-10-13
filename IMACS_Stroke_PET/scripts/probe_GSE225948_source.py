"""Auditable GEO GSE225948 source schema probe. No expression inference."""
import csv, gzip, io, json, re, hashlib
from pathlib import Path
import requests
from bs4 import BeautifulSoup

root=Path("IMACS_Stroke_PET/results/source_probes")
root.mkdir(parents=True,exist_ok=True)
ids=["GSM7060815","GSM7060816","GSM7060817","GSM7060818","GSM7060819","GSM7060820","GSM7060821","GSM7060822"]
rows=[]
session=requests.Session()
for gsm in ids:
    url=f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/"
    response=session.get(url,timeout=120)
    print(gsm,response.status_code,url,flush=True)
    response.raise_for_status()
    soup=BeautifulSoup(response.text,"html.parser")
    names=[a.get("href","") for a in soup.find_all("a")]
    names=[n for n in names if n.startswith(gsm+"_") and (n.endswith("_counts.csv.gz") or n.endswith("_metadata.csv.gz"))]
    matches={"counts":[n for n in names if n.endswith("_counts.csv.gz")],"metadata":[n for n in names if n.endswith("_metadata.csv.gz")]}
    for kind,candidates in matches.items():
        if len(candidates)!=1:
            raise ValueError(f"{gsm} expected one {kind} file, found {candidates}. Directory had: {names[:10]}")
        fname=candidates[0]
        if kind=="metadata":
            file_url=url+fname
            binary=session.get(file_url,timeout=180);binary.raise_for_status()
            raw=binary.content
            text=gzip.decompress(raw).decode("utf-8-sig")
            sample=pd.read_csv(io.StringIO(text),nrows=6) if False else None
            first=list(csv.reader(io.StringIO(text)))[:6]
            schema={"header":first[0],"first_rows":first[1:4],"total_csv_lines":len(text.splitlines())}
            (root/f"{gsm}_metadata_schema.json").write_text(json.dumps(schema,indent=2))
            rows.append({"gsm":gsm,"kind":kind,"url":file_url,"filename":fname,"bytes":len(raw),
                         "sha256":hashlib.sha256(raw).hexdigest(),"n_columns_firstline":len(first[0]),
                         "columns_firstline":";".join(first[0][:20])})
        else:
            file_url=url+fname
            rr=session.get(file_url,stream=True,timeout=240);rr.raise_for_status()
            gf=gzip.GzipFile(fileobj=rr.raw)
            header=gf.readline().decode("utf-8-sig").rstrip("\r\n")
            firstline=gf.readline().decode("utf-8-sig").rstrip("\r\n")
            rr.close()
            h=next(csv.reader([header]))
            first=next(csv.reader([firstline]))
            (root/f"{gsm}_counts_schema.json").write_text(json.dumps({"header_first20":h[:20],
                 "header_count":len(h),"first_gene_row_20":first[:20],"first_row_length":len(first)},indent=2))
            rows.append({"gsm":gsm,"kind":kind,"url":file_url,"filename":fname,
                         "bytes":None,"sha256":None,"n_columns_firstline":len(h),
                         "columns_firstline":";".join(h[:20])})
import pandas as pd
pd.DataFrame(rows).to_csv(root/"GSE225948_day2_input_source_manifest.csv",index=False)
print("GSE225948 schema probe completed",len(rows),"files.",flush=True)
