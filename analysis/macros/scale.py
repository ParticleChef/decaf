import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import processor 
import hist
from coffea.util import load, save

def scale(filename):

    hists = load(filename)

    ###
    # Rescaling MC histograms using the xsec weight
    ###

    scale={}
    for dataset in hists['sumw'].keys():
        scale[dataset]=hists['sumw'][dataset]
    print('Sumw extracted')

    for key in hists.keys():
        if key=='sumw': continue
        for dataset in hists[key].keys():
            if 'MET' in dataset or 'SingleElectron' in dataset or 'SinglePhoton' in dataset or 'EGamma' in dataset or 'BTagMu' in dataset: continue
            hists[key][dataset] *= 1/scale[dataset]
    print('Histograms scaled')


    ###
    # Defining 'process', to aggregate different samples into a single process
    ##
    sig_map = {}
    bkg_map = {}
    data_map = {}
    #bkg_map["VV"] = (["WW*","WZ*","ZZ*"],)
    #bkg_map["WW"] = ["WW"]
    #bkg_map["WZ"] = ["WZ"]
    #bkg_map["ZZ"] = ["ZZ"]
    bkg_map["QCD Multijet"] = ["QCD"]
    bkg_map["VV"] = ["_TuneCP5_13TeV-pythia8"]
    bkg_map["TT"] = ["TT"]
    bkg_map["Single Top"] = ["ST"]
    bkg_map[r"W ($\ell\nu$) + Jets"]  = ["WJets"]
    bkg_map[r"Z ($\ell\ell$) + Jets"] = ["DYJets"]
    bkg_map[r"Z ($\nu\nu$) + Jets"]   = ["JetsToNuNu"]
    bkg_map["W+HF"] = ["HF--WJets"]
    bkg_map["W+LF"] = ["LF--WJets"]
    bkg_map["W+jetsLO"] = ["lo--WJets"]
    bkg_map["W+jetsNNLO"] = ["nnlo--WJets"]
    bkg_map["Z+HF"] = ["HF--ZJetsToNuNu"]
    bkg_map["Z+LF"] = ["LF--ZJetsToNuNu"]
    bkg_map["Z+jetsLO"] = ["lo--ZJets"]
    bkg_map["Z+jetsNNLO"] = ["nnlo--ZJets"]
    bkg_map["G+HF"] = ["HF--GJets"]
    bkg_map["G+LF"] = ["LF--GJets"]
    bkg_map["G+jetsLO"] = ["lo--GJets"]
    bkg_map["G+jetsNNLO"] = ["nnlo--GJets"]
    data_map["MET"] = ["MET"]
    data_map["SingleElectron"] = ["SingleElectron"]
    data_map["SinglePhoton"] = ["SinglePhoton"]
    data_map["EGamma"] = ["EGamma"]
    data_map["BTagMu"] = ["BTagMu"]

    for signal in hists['sumw'].keys():
        if 'TPhiTo2Chi' not in signal: continue
        print(signal)
        sig_map[signal] = signal  ## signals
    print('Processes defined')
    
    ###
    # Storing signal and background histograms
    ###
    bkg_hists={}
    sig_hists={}
    data_hists={}
    for key in hists.keys():
        bkg_hists[key]={}
        sig_hists[key]={}
        data_hists[key]={}
        for process in bkg_map.keys():
            for dataset in hists[key].keys():
                if not any(d in dataset for d in bkg_map[process]): continue
                #print('Adding',dataset,'to',process,'for variable',key)
                try:
                    bkg_hists[key][process]+=hists[key][dataset]
                except:
                    bkg_hists[key][process]=hists[key][dataset]
        for process in data_map.keys():
            for dataset in hists[key].keys():
                if not any(d in dataset for d in data_map[process]): continue
                #print('Adding',dataset,'to',process,'for variable',key)
                try:
                    data_hists[key][process]+=hists[key][dataset]
                except:
                    data_hists[key][process]=hists[key][dataset]
        for process in sig_map.keys():
            for dataset in hists[key].keys():
                if not any(d in dataset for d in sig_map[process]): continue
                #print('Adding',dataset,'to',process,'for variable',key)
                try:
                    sig_hists[key][process]+=hists[key][dataset]
                except:
                    sig_hists[key][process]=hists[key][dataset]
        #for signal in sig_hists[key].keys():
        #    print('Scaling '+ signal +' by xsec '+str(xsec[signal]))
        #    sig_hists[key] *= xsec[str(signal)]
        
    print('Histograms grouped')

    return bkg_hists, sig_hists, data_hists

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--file', help='file', dest='file')
    (options, args) = parser.parse_args()

    bkg_hists, sig_hists, data_hists = scale(options.file)
    name = options.file

    hists={
        'bkg': bkg_hists,
        'sig': sig_hists,
        'data': data_hists
    }
    save(hists,name.replace('.merged','.scaled'))
