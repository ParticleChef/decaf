import json
import os
import glob
import gzip
import os, sys

year = '2024'
outdict = {}
      
### Data Loop
# search directory
if year == '2022pre':
    storage_dirlist = ['/data/data/2022C/', '/data/data/2022D/']
    mc_dir = '/data/mc/Run3Summer22NanoAODv12'
elif year == '2022post':
    storage_dirlist = ['/data/data/2022E/', '/data/data/2022F/', '/data/data/2022G/']
    mc_dir = '/data/mc/Run3Summer22EENanoAODv12'
elif year == '2023pre':
    storage_dirlist = ['/data/data/2023C/']
    mc_dir = '/data/mc/Run3Summer23NanoAODv12'
elif year == '2023post':
    storage_dirlist = ['/data/data/2023D/']
    mc_dir = '/data/mc/Run3Summer23BPixNanoAODv12'
elif year == '2024':
    storage_dirlist = ['/scratch/twkim/data/2024C/', '/scratch/twkim/data/2024D/', '/scratch/twkim/data/2024E/', '/scratch/twkim/data/2024F/', '/scratch/twkim/data/2024G/', '/scratch/twkim/data/2024H/', '/scratch/twkim/data/2024I/']
    data_dirlist = ['/data/data/2024C/', '/data/data/2024D/', '/data/data/2024E/', '/data/data/2024F/', '/data/data/2024G/', '/data/data/2024H/', '/data/data/2024I/']
    mc_dir = '/data/mc/RunIII2024Summer24NanoAODv15'
    #mc_dir = '/data/mc/stop_2024_nTuple'

for storage_dir in storage_dirlist:
    # get ls
    ls = os.listdir(storage_dir)
    for key in ls:
        ## EGamma pass
        if 'EGamma' in key:
            continue
        files = os.listdir(storage_dir+'/'+key)
        files = [storage_dir+'/'+key+'/'+f for f in files]
        if len(files) > 100:
            idx = 1
            for i in range(0, len(files), 100):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+100]
                outdict.update({new_key: {'files': new_data, 'xs': -1}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': -1}})

for data_dir in data_dirlist:
    # get ls
    ls = os.listdir(data_dir)
    for key in ls:
        ## JetMET pass
        if 'JetMET' in key:
            continue
        if key == 'EGamma0' or key == 'EGamma1':
            continue
        files = os.listdir(data_dir+'/'+key)
        files = [data_dir+'/'+key+'/'+f for f in files]
        if len(files) > 20:
            idx = 1
            for i in range(0, len(files), 20):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+20]
                outdict.update({new_key: {'files': new_data, 'xs': -1}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': -1}})

# MC data
mc_ls = os.listdir(mc_dir)
for key in mc_ls:
    if 'backup' in key:
        continue
    print('Processing MC key:', key)
    files = os.listdir(mc_dir+'/'+key)
    files = [mc_dir+'/'+key+'/'+f for f in files]
    ## TTto4Q
    if 'TTto4Q' in key:
        if len(files) > 40:
            idx = 1
            for i in range(0, len(files), 40):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+40]
                outdict.update({new_key: {'files': new_data, 'xs': 419.82}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 419.82}})
    ## TTtoLNu2Q
    elif 'TTtoLNu2Q' in key:
        if len(files) > 40:
            idx = 1
            for i in range(0, len(files), 40):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+40]
                outdict.update({new_key: {'files': new_data, 'xs': 405.75}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 405.75}})
    #TTto2L2Nu
    elif 'TTto2L2Nu' in key:
        if len(files) > 40:
            idx = 1
            for i in range(0, len(files), 40):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+40]
                outdict.update({new_key: {'files': new_data, 'xs': 98.04}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 98.04}})
    ## QCD_Bin-PT-15to20
    elif 'QCD_Bin-PT-15to20' in key:
        if len(files) > 100:
            idx = 1
            for i in range(0, len(files), 100):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+100]
                outdict.update({new_key: {'files': new_data, 'xs': 885700000.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 885700000.0}})
    ## QCD_Bin-PT-20to30
    elif 'QCD_Bin-PT-20to30' in key:
        if len(files) > 100:
            idx = 1
            for i in range(0, len(files), 100):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+100]
                outdict.update({new_key: {'files': new_data, 'xs': 415700000.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 415700000.0}})
    ## QCD_Bin-PT-30to50
    elif 'QCD_Bin-PT-30to50' in key:
        if len(files) > 100:
            idx = 1
            for i in range(0, len(files), 100):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+100]
                outdict.update({new_key: {'files': new_data, 'xs': 112300000.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 112300000.0}})
    ## QCD_Bin-PT-50to80
    elif 'QCD_Bin-PT-50to80' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 16730000.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 16730000.0}})
    ## QCD_Bin-PT-80to120
    elif 'QCD_Bin-PT-80to120' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 2506000.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 2506000.0}})
    ## QCD_Bin-PT-120to170
    elif 'QCD_Bin-PT-120to170' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 439800.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 439800.0}})
    ## QCD_Bin-PT-170to300
    elif 'QCD_Bin-PT-170to300' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 113300.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 113300.0}})
    ## QCD_Bin-PT-300to470
    elif 'QCD_Bin-PT-300to470' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 7581.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 7581.0}})
    ## QCD_Bin-PT-470to600
    elif 'QCD_Bin-PT-470to600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 623.3}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 623.3}})
    ## QCD_Bin-PT-600to800
    elif 'QCD_Bin-PT-600to800' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 178.7}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 178.7}})
    ## QCD_Bin-PT-800to1000
    elif 'QCD_Bin-PT-800to1000' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 30.62}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 30.62}})
    ## QCD_Bin-PT-1000to1500
    elif 'QCD_Bin-PT-1000to1500' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 9.306}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 9.306}})
    ## QCD_Bin-PT-1500to2000
    elif 'QCD_Bin-PT-1500to2000' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.5015}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.5015}})
    ## QCD_Bin-PT-2000to2500
    elif 'QCD_Bin-PT-2000to2500' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.04264}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.04264}})
    ## QCD_Bin-PT-2500to3000
    elif 'QCD_Bin-PT-2500to3000' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.004454}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.004454}})
    ## QCD_Bin-PT-3000
    elif 'QCD_Bin-PT-3000' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.0005539}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.0005539}})

    ## Zto2Nu-2Jets_Bin-1J-PTNuNu-40to100
    elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-40to100' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 885.7}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 885.7}})
    ## Zto2Nu-2Jets_Bin-1J-PTNuNu-100to200
    elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-100to200' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 78.91}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 78.91}})
    ## Zto2Nu-2Jets_Bin-1J-PTNuNu-200to400
    elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-200to400' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 5.406}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 5.406}})
    ## Zto2Nu-2Jets_Bin-1J-PTNuNu-400to600
    elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-400to600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.1693}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.1693}})
    ## Zto2Nu-2Jets_Bin-1J-PTNuNu-600
    elif 'Zto2Nu-2Jets_Bin-1J-PTNuNu-600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.01895}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.01895}})
    ## Zto2Nu-2Jets_Bin-2J-PTNuNu-40to100
    elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-40to100' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 328.6}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 328.6}})
    ## Zto2Nu-2Jets_Bin-2J-PTNuNu-100to200
    elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-100to200' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 100.2}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 100.2}})
    ## Zto2Nu-2Jets_Bin-2J-PTNuNu-200to400
    elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-200to400' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 13.76}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 13.76}})
    ## Zto2Nu-2Jets_Bin-2J-PTNuNu-400to600
    elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-400to600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.775}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.775}})
    ## Zto2Nu-2Jets_Bin-2J-PTNuNu-600
    elif 'Zto2Nu-2Jets_Bin-2J-PTNuNu-600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.1304}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.1304}})
    elif 'WtoLNu-2Jets_Bin-1J-PTLNu-40to100' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 4211.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 4211.0}})
    elif 'WtoLNu-2Jets_Bin-1J-PTLNu-100to200' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 342.3}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 342.3}})
    elif 'WtoLNu-2Jets_Bin-1J-PTLNu-200to400' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 21.84}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 21.84}})
    elif 'WtoLNu-2Jets_Bin-1J-PTLNu-400to600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.6845}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.6845}})
    elif 'WtoLNu-2Jets_Bin-1J-PTLNu-600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.07753}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.07753}})
    elif 'WtoLNu-2Jets_Bin-2J-PTLNu-40to100' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 1581.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 1581.0}})
    elif 'WtoLNu-2Jets_Bin-2J-PTLNu-100to200' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 411.1}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 411.1}})
    elif 'WtoLNu-2Jets_Bin-2J-PTLNu-200to400' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 53.59}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 53.59}})
    elif 'WtoLNu-2Jets_Bin-2J-PTLNu-400to600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 3.099}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 3.099}})
    elif 'WtoLNu-2Jets_Bin-2J-PTLNu-600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.5259}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.5259}})
            
    elif 'ST-TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 3.848}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 3.848}})
    elif 'ST-TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 16.338}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 16.338}})
    elif 'ST-TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 15.857}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 15.857}})
    elif 'ST-TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 3.841}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 3.841}})
    elif 'ST-TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 16.311}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 16.311}})
    elif 'ST-TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 15.831}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 15.831}})
    elif 'ST-TBbarQto2Q-t-channel-4FS' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 77.26}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 77.26}})
    elif 'ST-TBbarQtoLNu-t-channel-4FS' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 38.6}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 38.6}})
    elif 'ST-TbarBQto2Q-t-channel-4FS' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 46.73}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 46.73}})
    elif 'ST-TbarBQtoLNu-t-channel-4FS' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 23.34}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 23.34}})
    elif 'SMS-2Stop_Par-mStop-600' in key:
        if len(files) > 100:
            idx = 1
            for i in range(0, len(files), 100):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+100]
                outdict.update({new_key: {'files': new_data, 'xs': 10.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 10.0}})
    elif 'SMS-2Stop_Par-mStop-1000' in key:
        if len(files) > 100:
            idx = 1
            for i in range(0, len(files), 100):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+100]
                outdict.update({new_key: {'files': new_data, 'xs': 10.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 10.0}})
    elif 'SMS-2Stop_Par-mStop-1500' in key:
        if len(files) > 100:
            idx = 1
            for i in range(0, len(files), 100):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+100]
                outdict.update({new_key: {'files': new_data, 'xs': 10.0}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 10.0}})
    elif 'TBbartoLplusNuBbar-s-channel-4FS' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 2.278}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 2.278}})
    elif 'TbarBtoLminusNuB-s-channel-4FS' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 1.43}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 1.43}})
    elif 'VV-WW_TuneCP5_13p6TeV_pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 80.23}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 80.23}})
    elif 'VV-WZ_TuneCP5_13p6TeV_pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 29.1}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 29.1}})
    elif 'VV-ZZ_TuneCP5_13p6TeV_pythia8' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 12.75}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 12.75}})
    elif 'GJ_Bin-PTG-100to200' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 1391}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 1391}})
    elif 'GJ_Bin-PTG-200to400' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 88.24}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 88.24}})
    elif 'GJ_Bin-PTG-400to600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 3.77}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 3.77}})
    elif 'GJ_Bin-PTG-600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.576}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.576}})
    elif 'DYto2L-2Jets_Bin-MLL-50-PTLL-100' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 106.6}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 106.6}})
    elif 'DYto2L-2Jets_Bin-MLL-50-PTLL-200' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 10.83}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 10.83}})
    elif 'DYto2L-2Jets_Bin-MLL-50-PTLL-400' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.5943}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.5943}})
    elif 'DYto2L-2Jets_Bin-MLL-50-PTLL-600' in key:
        if len(files) > 50:
            idx = 1
            for i in range(0, len(files), 50):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+50]
                outdict.update({new_key: {'files': new_data, 'xs': 0.08108}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': 0.08108}})

# update json file
with open('KNU_'+ str(year) +'_v2.json', 'w') as f:
    json.dump(outdict, f, indent=4)
# recompress the json file
os.system('gzip -f KNU_'+ str(year) +'_v2.json')
