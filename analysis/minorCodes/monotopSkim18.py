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
def analisi(filelist,samplename, max_events = 200000, isData=True):
    flist = []
    for f in filelist:
        flist.append(f+":Events")
    
    exceptlist = ['GenMET', 'TkMET', 'ChsMET', 'RawMET','AK4PFPuppi_Jet', 'FatJet',
                    'TrigObj', 'SV', 'PPSLocalTrack', 'TrigObj', 'Proton', 'SubJet',
                    'GenIsolatedPhoton', 'GenDressedLepton']

    storelist = ['run', 'luminosityBlock', 'event', 'L1PreFiringWeight', 'fixedGridRhoFastjetAll', 'PV',
                    'nElectron', 'nAK15PFPuppi', 'nAK15PFPuppi', 'nJet', 'nMuon', 'nPhoton', 'nTau',
                    'Electron', 'AK15PFPuppi', 'AK15PFPuppi', 'Jet', 'MET', 'Muon', 'Photon', 'Tau',
                    'Flag' ]

    branches = []
    br = uproot.open(flist[0])
    for i in br.keys():
        if 'HLT' in i:
            if "HLT_Ele32_WPTight_Gsf" in i or "HLT_Photon200" in i or "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60" in i or "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight" in i:
                branches.append(i)
                print('store branch: ', i)
        elif i.split('_')[0] in storelist:
            branches.append(i)
            print('store branch: ', i)
        else:
            continue

    number_of_events = 0
    skimmed_events = 0
    all_arrays = None
    idx = 0

    for bb in branches:
        if 'n' in bb:
            print(bb)

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

        for bb in branches:
            if 'n' in bb:
                print('2 ', bb)
        while len(all_arrays) >= max_events:
            chunk = all_arrays[:max_events]
            all_arrays = all_arrays[max_events:]

        # make output directory
            outdir = outputdir + key + "/"
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            # save root file
            num = str(idx)
            with uproot.recreate(outdir+samplename+"_skimmed_"+num+".root") as f:
                #for branch in branches:
                #    data = ak.to_numpy(chunk[branch]) if isinstance(chunk[branch], ak.Array) else chunk[branch]
                #    f["Events"] = {branch: data}
                f["Events"] = {
                    #branch: chunk[branch].to_numpy() if isinstance(chunk[branch], ak.Array) else chunk[branch]  for branch in branches
                    branch: chunk[branch] for branch in branches
                }
#            with uproot.update(outdir+samplename+"_skimmed_"+num+".root") as ff:
#                if 'nElectron_pt' in ff["Events"]:
#                    del ff["Events/nElectron_pt"]
#            #    for branchname in uproot.open(outdir+samplename+"_skimmed_"+num+".root"):
#            #        if not branchname in branches:
#            #            delbranch = "Events/"+str(branchname)
#            #            del ff[delbranch]

            #remove_branch(outdir+samplename+"_skimmed_"+num+".root", "nElectron_pt")
            idx += 1
        for bb in branches:
            if 'n' in bb:
                print('2 ', bb)

    if len(all_arrays) > 0:
        outdir = outputdir + key + "/"
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        num = str(idx)
        with uproot.recreate(outdir+samplename+"_skimmed_"+num+".root") as f:
            #for branch in branches:
            #    data = ak.to_numpy(chunk[branch]) if isinstance(chunk[branch], ak.Array) else chunk[branch]
            #    f["Events"] = {branch: data}
            f["Events"] = {
                #branch: all_arrays[branch].to_numpy() if isinstance(all_arrays[branch], ak.Array) else all_arrays[branch]   for branch in branches
                branch: chunk[branch] for branch in branches
            }

#        with uproot.update(outdir+samplename+"_skimmed_"+num+".root") as ff:
#            if 'nElectron_pt' in ff["Events"]:
#                del ff["Events/nElectron_pt"]
#        #    for branchname in uproot.open(outdir+samplename+"_skimmed_"+num+".root"):
#        #        if not branchname in branches:
#        #            delbranch = "Events/"+str(branchname)
#        #            del ff[delbranch]
        #remove_branch(outdir+samplename+"_skimmed_"+num+".root", "nElectron_pt")

#def remove_branch(filepath, branch_name):
#    with uproot.open(filepath) as file:
#        events = file["Events"]
#        arrays = events.arrays()
#        new_arrays = {key: value for key, value in arrays.items() if key != branch_name}
#
#        with uproot.recreate(filepath) as new_file:
#            new_file["Events"] = arrays

metadata = "/home/jhong/monotop/CMSSW_11_3_4/src/decaf/analysis/metadata/onefile.json.gz"
#metadata = "/home/jhong/jieun_mtop_decaf/analysis/metadata/KNUv1_UL_ALL2018_v3.json"
outputdir = "/data/mTopnTuples/monotop2018/"
#keys = ['WJetsToLNu_Pt-600ToInf']
keys = ['WJetsToLNu_Pt-400To600']


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
