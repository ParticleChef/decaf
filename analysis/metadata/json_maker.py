import json
import os
import glob
import gzip
import os, sys

year = '2022pre'
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

for storage_dir in storage_dirlist:
    # get ls
    ls = os.listdir(storage_dir)
    for key in ls:
        files = os.listdir(storage_dir+'/'+key)
        files = [storage_dir+'/'+key+'/'+f for f in files]
        if len(files) > 10:
            idx = 1
            for i in range(0, len(files), 10):
                new_key = key + '____' + str(idx) + '_'
                new_data = files[i:i+10]
                outdict.update({new_key: {'files': new_data, 'xs': -1}})
                idx += 1
        else:
            outdict.update({key: {'files': files, 'xs': -1}})

# MC data
mc_ls = os.listdir(mc_dir)
for key in mc_ls:
    # skip QCD
    if 'QCD' in key or 'backup' in key:
        continue
    print('Processing MC key:', key)
    files = os.listdir(mc_dir+'/'+key)
    files = [mc_dir+'/'+key+'/'+f for f in files]
    if len(files) > 10:
        idx = 1
        for i in range(0, len(files), 10):
            new_key = key + '____' + str(idx) + '_'
            new_data = files[i:i+10]
            outdict.update({new_key: {'files': new_data, 'xs': -1}})
            idx += 1
    else:
        outdict.update({key: {'files': files, 'xs': -1}})

# update json file
with open('KNU_'+ str(year) +'_v1.json', 'w') as f:
    json.dump(outdict, f, indent=4)
# recompress the json file
os.system('gzip -f KNU_'+ str(year) +'_v1.json')


