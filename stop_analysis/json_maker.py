import json
import os
import glob

year = 2023
outdict = {}
      
### Data Loop
# search directory
storage_dirlist = glob.glob('/data/data/*'+ str(year) +'*/')
#['/data/data/2018A/', '/data/data/2018B/', '/data/data/2018C/', '/data/data/2018D/']
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

print(outdict)
# sorting by key

# update json file
with open('KNU_'+ str(year) +'_v1.json', 'w') as f:
    json.dump(outdict, f, indent=4)

