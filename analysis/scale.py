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
        #scale[dataset]=hists['sumw'][dataset]
        ### This scale is only for 2024 samples
        '''
        if 'TTto4Q' in dataset:
            scale[dataset] = 164126770000
        elif 'TTtoLNu2Q' in dataset:
            scale[dataset] = 161706980000
        #elif 'TTto2L2Nu' in dataset:
        #    scale[dataset] = 37820105000
        elif 'QCD_Bin-PT-15to20' in dataset:
            scale[dataset] = 99015416
        elif 'QCD_Bin-PT-20to30' in dataset:
            scale[dataset] = 93237180
        elif 'QCD_Bin-PT-30to50' in dataset:
            scale[dataset] = 99983304
        elif 'QCD_Bin-PT-50to80' in dataset:
            scale[dataset] = 97814970
        elif 'QCD_Bin-PT-80to120' in dataset:
            scale[dataset] = 99253510
        elif 'QCD_Bin-PT-120to170' in dataset:
            scale[dataset] = 99957400
        elif 'QCD_Bin-PT-170to300' in dataset:
            scale[dataset] = 99860960
        elif 'QCD_Bin-PT-300to470' in dataset:
            scale[dataset] = 75346750
        elif 'QCD_Bin-PT-470to600' in dataset:
            scale[dataset] = 77529080
        elif 'QCD_Bin-PT-600to800' in dataset:
            scale[dataset] = 76078900
        elif 'QCD_Bin-PT-800to1000' in dataset:
            scale[dataset] = 79492040
        elif 'QCD_Bin-PT-1000to1500' in dataset:
            scale[dataset] = 79876520
        elif 'QCD_Bin-PT-1500to2000' in dataset:
            scale[dataset] = 19997308
        elif 'QCD_Bin-PT-2000to2500' in dataset:
            scale[dataset] = 19311060
        elif 'QCD_Bin-PT-2500to3000' in dataset:
            scale[dataset] = 19988864
        elif 'QCD_Bin-PT-3000toInf' in dataset:
            scale[dataset] = 19347696
        elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-40to100' in dataset:
            scale[dataset] = 2113813300000
        elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-100to200' in dataset:
            scale[dataset] = 155959030000
        elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-200to400' in dataset:
            scale[dataset] = 879827600
        elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-400to600' in dataset:
            scale[dataset] = 26588450
        elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-600' in dataset:
            scale[dataset] = 3285069.5
        elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-40to100' in dataset:
            scale[dataset] = 1053125200000
        elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-100to200' in dataset:
            scale[dataset] = 173473170000
        elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-200to400' in dataset:
            scale[dataset] = 2170973200
        elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-400to600' in dataset:
            scale[dataset] = 29346958
        elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-600' in dataset:
            scale[dataset] = 3812532
        elif 'WtoLNu-2Jets_Bin-1J-PTLNu-40to100' in dataset:
            scale[dataset] = 20889200000000
        elif 'WtoLNu-2Jets_Bin-1J-PTLNu-100to200' in dataset:
            scale[dataset] = 1124590700000
        elif 'WtoLNu-2Jets_Bin-1J-PTLNu-200to400' in dataset:
            scale[dataset] = 13632560000
        elif 'WtoLNu-2Jets_Bin-1J-PTLNu-400to600' in dataset:
            scale[dataset] = 91151064
        elif 'WtoLNu-2Jets_Bin-1J-PTLNu-600' in dataset:
            scale[dataset] = 13635438
        elif 'WtoLNu-2Jets_Bin-2J-PTLNu-40to100' in dataset:
            scale[dataset] = 10847134000000
        elif 'WtoLNu-2Jets_Bin-2J-PTLNu-100to200' in dataset:
            scale[dataset] = 1192042600000
        elif 'WtoLNu-2Jets_Bin-2J-PTLNu-200to400' in dataset:
            scale[dataset] = 6849107000
        elif 'WtoLNu-2Jets_Bin-2J-PTLNu-400to600' in dataset:
            scale[dataset] = 105725900
        elif 'WtoLNu-2Jets_Bin-2J-PTLNu-600' in dataset:
            scale[dataset] = 14914732
        elif 'TWminusto2L2Nu' in dataset:
            scale[dataset] = 56972810
        elif 'TWminusto4Q' in dataset:
            scale[dataset] = 393607230
        elif 'TWminustoLNu2Q' in dataset:
            scale[dataset] = 1035079900
        elif 'TbarWplusto2L2Nu' in dataset:
            scale[dataset] = 57087280
        elif 'TbarWplusto4Q' in dataset:
            scale[dataset] = 394363400
        elif 'TbarWplustoLNu2Q' in dataset:
            scale[dataset] = 1059819970
        elif 'TBbarQto2Q-t-channel-4FS' in dataset:
            scale[dataset] = 6796076000
        elif 'TBbarQtoLNu-t-channel-4FS' in dataset:
            scale[dataset] = 1667160800
        elif 'TbarBQto2Q-t-channel-4FS' in dataset:
            scale[dataset] = 2025849000
        elif 'TbarBQtoLNu-t-channel-4FS' in dataset:
            scale[dataset] = 517007970
        else:
            scale[dataset]=hists['sumw'][dataset]
        '''
        scale[dataset]=hists['sumw'][dataset]
    #print('Sumw extracted')

    for key in hists.keys():
        if key=='sumw': continue
        #print(hists[key].keys())
        for dataset in hists[key].keys():
            print('Scaling',dataset,'for variable',key)
            if 'MET' in dataset or 'SingleElectron' in dataset or 'SinglePhoton' in dataset or 'EGamma' in dataset or 'BTagMu' in dataset or 'MuonEG' in dataset or 'Muon' in dataset or 'SingleMuon' in dataset: continue
            hists[key][dataset] *= 1/scale[dataset]
            print('Scaled',dataset,'for variable',key, 'by',scale[dataset])
    #print('Histograms scaled')


    ###
    # Defining 'process', to aggregate different samples into a single process
    ##

    sig_map = {}
    bkg_map = {}
    data_map = {}
    bkg_map['SMS-2Stop-Par-mStop-600'] = ["SMS-2Stop_Par-mStop-600"]
    bkg_map['SMS-2Stop-Par-mStop-1000'] = ["SMS-2Stop_Par-mStop-1000"]
    bkg_map['SMS-2Stop-Par-mStop-1500'] = ["SMS-2Stop_Par-mStop-1500"]
    bkg_map["QCD Multijet"] = ["QCD"]
    bkg_map["Z (inv)"] = ["Zto2Nu"]
    bkg_map["W (lnu)"] = ["WtoLNu"]
    bkg_map["Gamma + Jets"] = ["GJ"]
    bkg_map["VV"] = ["WW", "WZ", "ZZ"]
    #bkg_map["TT + V"] = ["TTZ","TTW"]
    bkg_map['TT'] = ["TTto", 'TTW', 'TTZ','TTTT','TTBB']
    bkg_map['DY'] = ['DY']
    #bkg_map["TT (AH)"] = ["TTto4Q"]
    #bkg_map["TT (SL)"] = ["TTtoLNu2Q"]
    #bkg_map["TT (DL)"] = ["TTto2L2Nu"]
    bkg_map["Single Top"] = [
        "TWminus",
        "TbarWplus",
        "TBbarQ",   # t-channel
        "TbarBQ",   # t-channel
        "TBbarto",  # s-channel
        "TbarBto",  # s-channel
    ]
    data_map["JetMET"] = ["JetMET"]
    data_map["Muon"] = ["Muon"]
    data_map["MuonEG"] = ["MuonEG"]
    data_map["SingleElectron"] = ["SingleElectron"]
    data_map["SinglePhoton"] = ["SinglePhoton"]
    data_map["EGamma"] = ["EGamma"]
    data_map["BTagMu"] = ["BTagMu"]
    #for signal in hists['sumw'].keys():
    #    if 'TTY2To2l2v2x' not in signal: continue
    #    print(signal)
    #    sig_map[signal] = signal  ## signals
    print('Processes defined')
    print(bkg_map.keys())
    
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
                if dataset != process: continue
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
