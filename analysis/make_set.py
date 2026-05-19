import json, gzip
meta="metadata/KNU_2024_v4.json.gz"
with gzip.open(meta, "rt") as f:
    d=json.load(f)
keys=sorted(d.keys())
with open("datasets.txt","w") as out:
    for k in keys:
        #if "Muon" in k:
            #out.write(k+"\n")
        out.write(k+"\n")
print("wrote datasets.txt with", len(keys), "datasets")