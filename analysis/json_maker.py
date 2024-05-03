import json
import os
from data.process import *

year = 2018

### MC Loop
# search directory
storage_dir = '/data/mc/UL2018NanoAODv9/'
outdict = {}
# get ls
ls = os.listdir(storage_dir)

for k, v in processes.items():
    if 'TTWJets' in k: continue
    if 'TTZToNuNu' in k: continue
    #if not 'TPhi' in k: continue
    print(k, v)
    if not k in ls: continue
    xs = v[2]
    files = os.listdir(storage_dir+'/'+str(k))
    files = [storage_dir+'/'+str(k)+'/'+f for f in files]
    if len(files) > 10:
        idx = 1
        for i in range(0,len(files),10):
            new_key = str(k)+'____'+str(idx)+'_'
            new_data = files[i:i+10]
            outdict.update({new_key:{'files':new_data,'xs':xs}})
            idx += 1
    else:
        outdict.update({k:{'files':files,'xs':xs}})


            
### Data Loop
# search directory
storage_dirlist = ['/data/data/2018A/', '/data/data/2018B/', '/data/data/2018C/', '/data/data/2018D/']
for storage_dir in storage_dirlist:
    # get ls
    ls = os.listdir(storage_dir)
    for key in ls:
        if 'MuonEG' in key: continue
        if 'DoubleMuon' in key: continue
        if 'SingleMuon' in key: continue
        # get file list with absolute path
        files = os.listdir(storage_dir+'/'+key)
        files = [storage_dir+'/'+key+'/'+f for f in files]
        if 'DoubleMuon' in key:
            if len(files) > 3:
                idx = 1
                for i in range(0, len(files), 3):
                    new_key = key + '____' + str(idx) + '_'
                    new_data = files[i:i+3]
                    outdict.update({new_key: {'files': new_data, 'xs': -1}})
                    idx += 1
            else:
                outdict.update({key: {'files': files, 'xs': -1}})
        elif 'MuonEG' in key:
            if len(files) > 3:
                idx = 1
                for i in range(0, len(files), 3):
                    new_key = key + '____' + str(idx) + '_'
                    new_data = files[i:i+3]
                    outdict.update({new_key: {'files': new_data, 'xs': -1}})
                    idx += 1
            else:
                outdict.update({key: {'files': files, 'xs': -1}})
        else:
            if len(files) > 10:
                idx = 1
                for i in range(0, len(files), 10):
                    new_key = key + '____' + str(idx) + '_'
                    new_data = files[i:i+10]
                    outdict.update({new_key: {'files': new_data, 'xs': -1}})
                    idx += 1
            else:
                outdict.update({key: {'files': files, 'xs': -1}})

print(outdict)
# sorting by key

# update json file
with open('KNUv1_UL_ALL2018_v3.json', 'w') as f:
    json.dump(outdict, f, indent=4)

