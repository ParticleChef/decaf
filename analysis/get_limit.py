import os, sys
import numpy as np
import matplotlib.pyplot as plt
import mplhep
import pandas as pd

# get all files in the directory
path = "datacards/hadmonotop2018_kps_edition/"
files = os.listdir("datacards/hadmonotop2018_kps_edition")

# solt files in order of MPhiXXX
files.sort(key=lambda x: int(x.split("_")[2].split("MPhi")[1].split(".")[0]))
df = pd.DataFrame(columns=["Mediator", "Dark Matter", "p2", "p1", "exp", "m1", "m2"])
for file in files:
    # get the mass point
    med = int(file.split("_")[2].split("MPhi")[1].split(".")[0])
    dm = int(file.split("_")[3].split("MChi")[1].split(".")[0])
    if dm*2 > med:
        print('skip...')
        continue
    print(f"Mediator mass: {med} GeV, Dark Matter mass: {dm} GeV")
    print(file)
    # run combine and get the limit
    os.system(f"combine -M AsymptoticLimits {path+file} --run blind > limit.txt")
    # get the limit

    with open("limit.txt") as f:
        lines = f.readlines()
        for idx,line in enumerate(lines):
            if "Expected" in line:
                try:
                    p2 = float(line.split(" ")[-1])
                    p1 = float(lines[idx+1].split(" ")[-1])
                    exp = float(lines[idx+2].split(" ")[-1])
                    m1 = float(lines[idx+3].split(" ")[-1])
                    m2 = float(lines[idx+4].split(" ")[-1])
                    break
                except:
                    continue
            
        # if no expected limit is found, skip the mass point
        try:
            df = df.append({"Mediator": med, "Dark Matter": dm, "p2": p2, "p1": p1, "exp": exp, "m1": m1, "m2": m2}, ignore_index=True)

            df.to_csv("limit.csv", index=False)
        except:
            continue
                

