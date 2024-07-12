import uproot
import numpy as np
import os
import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import multiprocessing as mp
import pandas as pd
import json
import time
import gzip

start_time = time.time()

### For 2017 dataset, same as 2018
# uproot iterate
def analisi(filelist,samplename,isData=True):
    idx = 0
    flist = []
    for f in filelist:
        flist.append(f+":Events")
    
    exceptlist = ['GenMET', 'TkMET', 'ChsMET', 'RawMET','AK4PFPuppi_Jet', 'FatJet',
                    'TrigObj', 'SV', 'PPSLocalTrack', 'TrigObj', 'Proton', 'SubJet',
                    'GenIsolatedPhoton', 'GenDressedLepton']
    branches = []
    br = uproot.open(flist[0])
    for i in br.keys():
        if 'HLT' in i:
            if "HLT_Ele32_WPTight_Gsf" in i or "HLT_Photon200" in i or "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60" in i or "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight" in i:
                branches.append(i)
        elif i.split('_')[0] in exceptlist:
            pass
        else:
            branches.append(i)
    print(branches)

    number_of_events = 0
    skimmed_events = 0
    all_arrays = None

    for arrays in uproot.iterate(flist,branches,entrysteps=30000):  ## entrysteps max is 40462 (around 20% cpu used hep)
        print(arrays)                                               ## Bigger number (ex. 10000000) automatically reduced
        number_of_events += len(arrays)
        print("processing events: ",number_of_events)
        
        ## MET filters
        goodVertices = arrays["Flag_goodVertices"]
        globalSuperTightHalo2016Filter = arrays["Flag_globalSuperTightHalo2016Filter"]
        HBHENoiseFilter = arrays["Flag_HBHENoiseFilter"]
        HBHENoiseIsoFilter = arrays["Flag_HBHENoiseIsoFilter"]
        EcalDeadCellTriggerPrimitiveFilter = arrays["Flag_EcalDeadCellTriggerPrimitiveFilter"]
        BadPFMuonFilter = arrays["Flag_BadPFMuonFilter"]
        BadPFMuonDzFilter = arrays["Flag_BadPFMuonDzFilter"]
        ecalBadCalibFilter = arrays["Flag_ecalBadCalibFilter"]
        met_filter = goodVertices & globalSuperTightHalo2016Filter & HBHENoiseFilter & HBHENoiseIsoFilter & EcalDeadCellTriggerPrimitiveFilter & BadPFMuonFilter & BadPFMuonDzFilter & ecalBadCalibFilter
        arrays = arrays[met_filter]

        # Trigger
        HLT_Ele32_WPTight_Gsf = arrays["HLT_Ele32_WPTight_Gsf"]
        HLT_Photon200 = arrays["HLT_Photon200"]
        HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = arrays["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"]
        HLT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = arrays["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"]
        trigger = HLT_Ele32_WPTight_Gsf | HLT_Photon200 | HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 | HLT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
        arrays = arrays[trigger]

        skimmed_events += len(arrays)
        print("skimmed events: ",skimmed_events)
        print("Skimming efficiency: ",(skimmed_events/number_of_events)*100,"%")

        if all_arrays is None:
            all_arrays = arrays
        else:
            all_arrays = ak.concatenate([all_arrays, arrays])

        # make output directory
    outdir = outputdir + key + "/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # save root file
    num = str(idx)
    idx += 1
    with uproot.create(outdir+samplename+"_skimmed_"+num+".root", compression=None) as f:
        outtuple = {}
        f["Events"] = {
            branch: ak.Array(all_arrays[branch]) for branch in branches
        }

    with uproot.open(outdir+samplename+"_skimmed_"+num+".root") as f:
        tr = f["Events"]
        print(tr.keys())
        #exit()


metadata = "/home/jhong/monotop/CMSSW_11_3_4/src/decaf/analysis/metadata/onefile.json.gz"
#metadata = "/home/jhong/jieun_mtop_decaf/analysis/metadata/KNUv1_UL_ALL2018_v3.json"
outputdir = "/data/mTopnTuples/monotop2018/"
keys = ['WJetsToLNu_Pt-600ToInf']


with gzip.open(metadata) as fin:
    samplefiles = json.load(fin)  

for dataset, info in samplefiles.items():
    for key in keys:
        if key in dataset:
            print('Processing:',dataset)
            files = info['files']
            samplename = str(dataset)
            print(samplename)
            analisi(filelist=files,samplename=samplename)

print('Job done! Jieun and Taiwoo are HAPPY!')
print('Total time: ' + str(time.time()-start_time) + ' seconds')
