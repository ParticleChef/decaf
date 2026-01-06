import json
import gzip
import os
from data.process import *

#####
## 10 root files are in 1 job
## The normal and 'ext' version of nanoaod are not seperated.
##
## For private 2022, 2023 sample (monotop)
## data: /data/data/privatedata and mc: /data/data/privatemc.
#####

year = 2022
name_key = '_private_v1'
json_name = str(year) + name_key # It will be json file name

nodata = False ## If 'False', data is also in list.

### MC Loop
# search directory
storage_dir = '/data/mc/privatemc/' + str(year)
print(storage_dir)
outdict = {}

# get ls
ls = os.listdir(storage_dir)

for k, v in processes.items():
    if 'Mphi' in k: continue ## Use if or if not to use specific dataset or exclude.
    #print(k, v)
    if not k in ls: continue
    xs = v[2]
    target_dir = os.listdir(storage_dir+'/'+str(k))
    if '.root' in target_dir[0]:
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

    else:
        ## Search root files in Lowest directory
        files = []
        for version_dir in target_dir:
            for dirpath, dirnames, filenames in os.walk(os.path.join(storage_dir, str(k), str(version_dir))):
                for f in filenames:
                    if f.endswith(".root"):
                        files.append(os.path.join(dirpath, f))

        if len(files) > 10:
            idx = 1
            for i in range(0, len(files), 10):
                new_key = str(k)+'____'+str(idx)+'_'
                new_data = files[i:i+10]
                outdict.update({new_key:{'files':new_data,'xs':xs}})
                idx += 1
        else:
            outdict.update({k:{'files':files,'xs':xs}})


            
### Data Loop
# search directory
storage_dirdata = {
    '2022': '/data/data/privatedata/2022'
}
for k, v in processes.items():
    if nodata: break
    if 'Mphi' in k: continue
    if not v[2] == 1: continue
    print(k, v)
    target_dir = os.listdir(storage_dirdata[str(year)]+'/'+str(k))
    if '.root' in target_dir[0]:
        files = os.listdir(storage_dirdata[str(year)]+'/'+str(k))
        files = [storage_dirdata[str(year)]+'/'+str(k)+'/'+f for f in files]
        if len(files) > 10:
            idx = 1
            for i in range(0,len(files),10):
                new_key = str(k)+'____'+str(idx)+'_'
                new_data = files[i:i+10]
                outdict.update({new_key:{'files':new_data,'xs':-1}})
                idx += 1
        else:
            outdict.update({k:{'files':files,'xs':-1}})

    else:
        ## Search root files in Lowest directory
        files = []
        for version_dir in target_dir:
            for dirpath, dirnames, filenames in os.walk(os.path.join(storage_dirdata[str(year)], str(k), str(version_dir))):
                for f in filenames:
                    if f.endswith(".root"):
                        files.append(os.path.join(dirpath, f))

        if len(files) > 10:
            idx = 1
            for i in range(0,len(files),10):
                new_key = str(k)+'____'+str(idx)+'_'
                new_data = files[i:i+10]
                outdict.update({new_key:{'files':new_data,'xs':-1}})
                idx += 1
        else:
            outdict.update({k:{'files':files,'xs':-1}})


#for storage_dir in storage_dirlist[str(year)]:
#    if nodata: break
#    # get ls
#    ls = os.listdir(storage_dir)
#    for key in ls:
#        print('key: ', key)
#        if 'MuonEG' in key: continue
#        if 'DoubleMuon' in key: continue
#        if 'SingleMuon' in key: continue
#        # get file list with absolute path
#        target_dir = os.listdir(storage_dir+'/'+str(key))
#        files = os.listdir(storage_dir+'/'+key)
#        files = [storage_dir+'/'+key+'/'+f for f in files]
#        if len(files) > 10:
#            idx = 1
#            for i in range(0, len(files), 10):
#                new_key = key + '____' + str(idx) + '_'
#                new_data = files[i:i+10]
#                outdict.update({new_key: {'files': new_data, 'xs': -1}})
#                idx += 1
#            else:
#                outdict.update({key: {'files': files, 'xs': -1}})

#print(outdict)
# sorting by key

with open(f'{json_name}.json', 'w') as f:
    json.dump(outdict, f, indent=4)
