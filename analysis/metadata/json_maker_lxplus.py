import json
import os
import glob
import gzip
import os, sys

def crawlDataset(dataset):
    os.system('dasgoclient --query="file dataset=%s" > files_%s.txt' % (dataset, dataset.split('/')[1]+dataset.split('/')[2].split('-')[0]))

year = '2024'
outdict = {}

datasets = [
    '/JetMET0/Run2024C-MINIv6NANOv15-v1/NANOAOD',
    '/JetMET0/Run2024D-MINIv6NANOv15-v1/NANOAOD',
    '/JetMET0/Run2024E-MINIv6NANOv15-v1/NANOAOD',
    '/JetMET0/Run2024F-MINIv6NANOv15-v2/NANOAOD',
    '/JetMET0/Run2024G-MINIv6NANOv15-v2/NANOAOD',
    '/JetMET0/Run2024H-MINIv6NANOv15-v2/NANOAOD',
    '/JetMET0/Run2024I-MINIv6NANOv15-v2/NANOAOD',
    '/JetMET0/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD',
    '/JetMET1/Run2024C-MINIv6NANOv15-v1/NANOAOD',
    '/JetMET1/Run2024D-MINIv6NANOv15-v1/NANOAOD',
    '/JetMET1/Run2024E-MINIv6NANOv15-v1/NANOAOD',
    '/JetMET1/Run2024F-MINIv6NANOv15-v2/NANOAOD',
    '/JetMET1/Run2024G-MINIv6NANOv15-v2/NANOAOD',
    '/JetMET1/Run2024H-MINIv6NANOv15-v2/NANOAOD',
    '/JetMET1/Run2024I-MINIv6NANOv15-v1/NANOAOD',
    '/JetMET1/Run2024I-MINIv6NANOv15_v2-v2/NANOAOD',
    '/EGamma0/Run2024C-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma0/Run2024D-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma0/Run2024E-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma0/Run2024F-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma0/Run2024G-MINIv6NANOv15-v2/NANOAOD',
    '/EGamma0/Run2024H-MINIv6NANOv15-v2/NANOAOD',
    '/EGamma0/Run2024I-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma0/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD',
    '/EGamma1/Run2024C-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma1/Run2024D-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma1/Run2024E-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma1/Run2024F-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma1/Run2024G-MINIv6NANOv15-v2/NANOAOD',
    '/EGamma1/Run2024H-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma1/Run2024I-MINIv6NANOv15-v1/NANOAOD',
    '/EGamma1/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD',
    '/Muon0/Run2024C-MINIv6NANOv15-v1/NANOAOD',
    '/Muon0/Run2024D-MINIv6NANOv15-v1/NANOAOD',
    '/Muon0/Run2024E-MINIv6NANOv15-v1/NANOAOD',
    '/Muon0/Run2024F-MINIv6NANOv15-v1/NANOAOD',
    '/Muon0/Run2024G-MINIv6NANOv15-v1/NANOAOD',
    '/Muon0/Run2024H-MINIv6NANOv15-v1/NANOAOD',
    '/Muon0/Run2024I-MINIv6NANOv15-v1/NANOAOD',
    '/Muon0/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD',
    '/Muon1/Run2024C-MINIv6NANOv15-v1/NANOAOD',
    '/Muon1/Run2024D-MINIv6NANOv15-v1/NANOAOD',
    '/Muon1/Run2024E-MINIv6NANOv15-v1/NANOAOD',
    '/Muon1/Run2024F-MINIv6NANOv15-v1/NANOAOD',
    '/Muon1/Run2024G-MINIv6NANOv15-v2/NANOAOD',
    '/Muon1/Run2024H-MINIv6NANOv15-v2/NANOAOD',
    '/Muon1/Run2024I-MINIv6NANOv15-v1/NANOAOD',
    '/Muon1/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD',
]

mcsets = {
    '/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 419.82,
    '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 405.75,
    '/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 98.04,
    '/QCD_Bin-PT-15to20_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 885700000,
    '/QCD_Bin-PT-20to30_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 415700000,
    '/QCD_Bin-PT-30to50_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 112300000,
    '/QCD_Bin-PT-50to80_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 16730000,
    '/QCD_Bin-PT-80to120_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 2506000,
    '/QCD_Bin-PT-120to170_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 439800,
    '/QCD_Bin-PT-170to300_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 113300,
    '/QCD_Bin-PT-300to470_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 7581,
    '/QCD_Bin-PT-470to600_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 623.3,
    '/QCD_Bin-PT-600to800_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 178.7,
    '/QCD_Bin-PT-800to1000_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 30.62,
    '/QCD_Bin-PT-1000to1500_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 9.306,
    '/QCD_Bin-PT-1500to2000_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.5015,
    '/QCD_Bin-PT-2000to2500_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.04264,
    '/QCD_Bin-PT-2500to3000_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.004454,
    '/QCD_Bin-PT-3000_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.0005539,
    #'/Zto2Nu-2Jets_Bin-PTNuNu-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 200.6,
    #'/Zto2Nu-2Jets_Bin-PTNuNu-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v1/NANOAODSIM': 20.41,
    #'/Zto2Nu-2Jets_Bin-PTNuNu-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 1.108,
    #'/Zto2Nu-2Jets_Bin-PTNuNu-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 0.1485,
    #'/WtoLNu-2Jets_Bin-PTLNu-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 863.6,
    #'/WtoLNu-2Jets_Bin-PTLNu-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v1/NANOAODSIM': 81.81,
    #'/WtoLNu-2Jets_Bin-PTLNu-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 4.453,
    #'/WtoLNu-2Jets_Bin-PTLNu-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v1/NANOAODSIM': 0.6027,
    '/Zto2Nu-2Jets_Bin-1J-PTNuNu-40to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 885.7,
    '/Zto2Nu-2Jets_Bin-1J-PTNuNu-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 78.91,
    '/Zto2Nu-2Jets_Bin-1J-PTNuNu-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 5.406,
    '/Zto2Nu-2Jets_Bin-1J-PTNuNu-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.1693,
    '/Zto2Nu-2Jets_Bin-1J-PTNuNu-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.01895,
    '/Zto2Nu-2Jets_Bin-2J-PTNuNu-40to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 328.6,
    '/Zto2Nu-2Jets_Bin-2J-PTNuNu-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 100.2,
    '/Zto2Nu-2Jets_Bin-2J-PTNuNu-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 13.76,
    '/Zto2Nu-2Jets_Bin-2J-PTNuNu-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.775,
    '/Zto2Nu-2Jets_Bin-2J-PTNuNu-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.1304,
    '/WtoLNu-2Jets_Bin-1J-PTLNu-40to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 4211,
    '/WtoLNu-2Jets_Bin-1J-PTLNu-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 342.3,
    '/WtoLNu-2Jets_Bin-1J-PTLNu-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 21.84,
    '/WtoLNu-2Jets_Bin-1J-PTLNu-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.6845,
    '/WtoLNu-2Jets_Bin-1J-PTLNu-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.07753,
    '/WtoLNu-2Jets_Bin-2J-PTLNu-40to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 1581,
    '/WtoLNu-2Jets_Bin-2J-PTLNu-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 411.1,
    '/WtoLNu-2Jets_Bin-2J-PTLNu-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 53.59,
    '/WtoLNu-2Jets_Bin-2J-PTLNu-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 3.099,
    '/WtoLNu-2Jets_Bin-2J-PTLNu-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.5259,
    '/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 3.841,
    '/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 16.311,
    '/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 15.831,
    '/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 3.848,
    '/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 16.338,
    '/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 15.857,
    '/TBbarQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 77.26,
    '/TBbarQtoLNu-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 38.6,
    '/TbarBQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 46.73,
    '/TbarBQtoLNu-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 23.34,
    '/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 2.278,
    '/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 1.43,
    '/GJ_Bin-PTG-100to200_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 1391,
    '/GJ_Bin-PTG-200to400_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 88.24,
    '/GJ_Bin-PTG-400to600_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 3.77,
    '/GJ_Bin-PTG-600_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.576,
    '/WW_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 80.23,
    '/WZ_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 29.1,
    '/ZZ_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 12.75,
    '/DYto2L-2Jets_Bin-MLL-50-PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 106.6,
    '/DYto2L-2Jets_Bin-MLL-50-PTLL-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v1/NANOAODSIM': 10.83,
    '/DYto2L-2Jets_Bin-MLL-50-PTLL-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.5943,
    '/DYto2L-2Jets_Bin-MLL-50-PTLL-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM': 0.08108,
    '/TTTT_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.009652,
    '/TTW-WtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.4678,
    '/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 0.6426,
    '/SMS-2Stop_Par-mStop-1000_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 9.123e-03,
    '/SMS-2Stop_Par-mStop-1500_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 3.912e-04,
    '/SMS-2Stop_Par-mStop-600_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM': 2.560e-01,
}

### Data Loop
# crawl datasets
for dataset in datasets:
    print('Crawling dataset: %s' % dataset)
    crawlDataset(dataset)
    # read file
    with open('files_'+dataset.split('/')[1]+dataset.split('/')[2].split('-')[0]+'.txt', 'r') as f:
        files = f.read().splitlines()
        if len(files) == 0:
            print('No files found for dataset: %s' % dataset)
            continue
        # all files must have 'root://cms-xrd-global.cern.ch/' prefix
        #files = ['root://eoscms.cern.ch//eos/cms/' + file for file in files]
        files = ['root://cms-xrd-global.cern.ch/' + file for file in files]
    if len(files) > 10:
        idx = 1
        for i in range(0, len(files), 10):
            new_key =  dataset.split('/')[1] + '-' + dataset.split('/')[2] + '____' + str(idx) + '_'
            new_data = files[i:i+10]
            outdict.update({new_key: {'files': new_data, 'xs': -1}})
            idx += 1
    else:
        new_key = dataset.split('/')[1] + '-' + dataset.split('/')[2]
        outdict.update({new_key: {'files': files, 'xs': -1}})
    # remove txt file
    os.remove('files_'+dataset.split('/')[1]+dataset.split('/')[2].split('-')[0]+'.txt')

### MC Loop
# crawl datasets
for dataset, xs in mcsets.items():
    print('Crawling dataset: %s' % dataset)
    crawlDataset(dataset)
    # read file
    with open('files_'+dataset.split('/')[1]+dataset.split('/')[2].split('-')[0]+'.txt', 'r') as f:
        files = f.read().splitlines()
        if len(files) == 0:
            print('No files found for dataset: %s' % dataset)
            continue
        # all files must have 'root://cms-xrd-global.cern.ch/' prefix
        #files = ['root://eoscms.cern.ch//eos/cms/' + file for file in files]
        files = ['root://cms-xrd-global.cern.ch/' + file for file in files]
    if len(files) > 10:
        idx = 1
        for i in range(0, len(files), 10):
            new_key =  dataset.split('/')[1] + '-' + dataset.split('/')[2] + '____' + str(idx) + '_'
            new_data = files[i:i+10]
            outdict.update({new_key: {'files': new_data, 'xs': xs}})
            idx += 1
    else:
        new_key = dataset.split('/')[1] + '-' + dataset.split('/')[2]
        outdict.update({new_key: {'files': files, 'xs': xs}})
    # remove txt file
    os.remove('files_'+dataset.split('/')[1]+dataset.split('/')[2].split('-')[0]+'.txt')

# update json file
with open('KNU_'+ str(year) +'_v4.json', 'w') as f:
    json.dump(outdict, f, indent=4)
# recompress the json file
os.system('gzip -f KNU_'+ str(year) +'_v4.json')
