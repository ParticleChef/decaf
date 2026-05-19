#!/usr/bin/env python3

import os
import sys

if len(sys.argv) != 3:
    print("Usage: check_futures.py <datasets.txt> <futures_dir>")
    sys.exit(1)

datasets_file = sys.argv[1]
futures_dir   = sys.argv[2]

# --- load datasets ---
with open(datasets_file) as f:
    datasets = [line.strip() for line in f if line.strip()]

# --- list futures files ---
futures_files = [
    f for f in os.listdir(futures_dir)
    if f.endswith(".futures")
]

futures_set = set(futures_files)

missing = []
found   = []
extra   = []

for ds in datasets:
    fname = f"{ds}.futures"
    if fname in futures_set:
        found.append(ds)
    else:
        missing.append(ds)

# extras: futures that are not in datasets.txt
dataset_set = set(datasets)
for f in futures_files:
    dsname = f.replace(".futures", "")
    if dsname not in dataset_set:
        extra.append(f)

# --- report ---
print(f"\nTotal datasets      : {len(datasets)}")
print(f"Total futures files : {len(futures_files)}")
print(f"Found               : {len(found)}")
print(f"Missing             : {len(missing)}")
print(f"Extra               : {len(extra)}")

with open("missing.txt", "w") as f:
    if missing:
        print("\n=== Missing datasets ===")
        for m in missing:
            print(m)
            f.write(m + "\n")


if extra:
    print("\n=== Extra futures (not in dataset list) ===")
    for e in extra:
        print(e)