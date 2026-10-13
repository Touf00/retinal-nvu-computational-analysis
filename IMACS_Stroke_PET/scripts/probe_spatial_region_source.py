"""Public-data provenance/availability audit for spatial region labels.
No credentials or protected download APIs. This probes metadata only.
"""
import json, os, sys, pathlib, requests
from datetime import datetime,timezone
root=pathlib.Path(os.environ.get("IMACS_SOURCE_PROBE_DIR","IMACS_Stroke_PET/results/source_probes"))
root.mkdir(parents=True,exist_ok=True)
endpoints=[
 "https://api.data.mendeley.com/datasets/publics/gnb2dsjms2/files?version=1",
 "https://api.data.mendeley.com/datasets/gnb2dsjms2/files?version=1",
 "https://api.data.mendeley.com/datasets/gnb2dsjms2?version=1",
 "https://api.mendeley.com/datasets/gnb2dsjms2?version=1",
]
out=[]
s=requests.Session()
s.headers.update({"User-Agent":"IMACS-Stroke-PET-Audit/1.0 research-contact","Accept":"application/json, application/vnd.mendeley-public-dataset.1+json"})
for url in endpoints:
 row={"url":url,"checked_utc":datetime.now(timezone.utc).isoformat()}
 try:
  r=s.get(url,timeout=45)
  row.update({"status":r.status_code,"content_type":r.headers.get("content-type",""),"bytes":len(r.content)})
  if r.ok:
   try:
    obj=r.json()
    if isinstance(obj,list):
     row["file_items"]=[{"filename":x.get("filename"),"id":x.get("id"),
                         "size":x.get("size"),"content_size":(x.get("content_details") or {}).get("size"),
                         "sha256":(x.get("content_details") or {}).get("sha256_hash"),
                         "has_download_url":bool((x.get("content_details") or {}).get("download_url"))}
                        for x in obj]
    elif isinstance(obj,dict):
     row["json_keys"]=list(obj)[:35]
     row["file_items"]=[{"filename":x.get("filename"),"id":x.get("id"),
                         "size":x.get("size")}
                      for x in obj.get("files",[]) if isinstance(x,dict)]
    else:row["response_class"]=type(obj).__name__
   except ValueError:row["response_preview"]=r.text[:250]
  else:row["response_preview"]=r.text[:160]
 except requests.RequestException as e:
  row["error_type"]=type(e).__name__
 out.append(row)
 print(url,row.get("status"),[(f.get("filename"),f.get("size")) for f in row.get("file_items",[])])
(root/"09_mendeley_public_region_source_probe.json").write_text(json.dumps(out,indent=2))
found=[x for x in out if x.get("file_items")]
if found:
 print("SOURCE INVENTORY: SUCCESS. Region annotations require downloading processed Seurat or obtaining smaller original author metadata.")
else:
 print("SOURCE INVENTORY: BLOCKED/UNAVAILABLE. Do not make lesion-location assertions.")
