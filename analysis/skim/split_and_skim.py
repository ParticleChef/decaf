import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import time
import json
import gzip

############
## Read json file
############
metadata = 'KNUv1_UL_ALL2018_v4'
with gzip.open('../metadata/'+metadata+'.json.gz') as fin:
    samplefiles = json.load(fin)

###########
## Make txt files for mTopSkim.C
###########
txtlist_dir = 'txtfile_list'
os.makedirs(txtlist_dir, exist_ok=True)

datasetlist = []
for dataset, values in samplefiles.items():
    if not 'ZZ_' in  dataset and not 'WZ_' in dataset: 
        continue
    onedataset = dataset.split('____')[0]
    if not onedataset in datasetlist:
        print("The first read of ", onedataset)
        datasetlist.append(onedataset)
        with open(f'{txtlist_dir}/{onedataset}.txt', 'w') as split_txt:
            for onefile in values['files']:
                #print("values['files']:", onefile)
                split_txt.write(onefile+'\n')
                #break
        split_txt.close()
    else:
        print(onedataset, ' is already created.')
        with open(f'{txtlist_dir}/{onedataset}.txt', 'a') as split_txt:
            for onefile in values['files']:
                split_txt.write(onefile+'\n')
        split_txt.close()

###########
## run mTopSkim.C ? 
###########
start_time = time.time()
for dataname in datasetlist:
    if not 'ZZ' in  dataname and not 'WZ' in dataname:
        continue
    if 'EGamma' in dataname or 'MET' in dataname:
        filetype = 'DATA'
    else:
        filetype = 'MC'
    print('Now ', dataname, ' starts skimming...')
    os.system(f'cp underconstruct_mTopSkim.C {dataname}_mTopSkim.C')
    os.system(f'sed -i "1 s/splitted_file_name/{dataname}/"  ./{dataname}_mTopSkim.C')
    os.system(f'sed -i "3 s/splitted_file_name/{txtlist_dir}\/{dataname}.txt/"  ./{dataname}_mTopSkim.C')
    os.system(f'sed -i "53 s/splitted_file_name/{dataname}/"  ./{dataname}_mTopSkim.C')
    os.system(f'sed -i "53 s/FILETYPE/{filetype}/"  ./{dataname}_mTopSkim.C')

    os.system(f'nohup root -l {dataname}_mTopSkim.C') ##?
    time.sleep(1)


print('Total time: ', str(time.time()-start_time), ' seconds')

a = '''

def run_macro(input_file, output_file, macro_file):
#    for i, input_file in enumerate(input_file_list):
    cmd = ['root', '-l', '-b', '-q', f'{macro_file}("{input_file}", "{output_file}", "Events")']
    subprocess.run(cmd)

if __name__ == "__main__":
    input_file = 'inputrootMichael.txt'   ## 
    lines_per_file = 10
    macro_file = 'mTopSkim.C'    ## 

    output_dir = '/10T/scratch/jhong/skimmingTest'
    os.makedirs(output_dir, exist_ok=True)

    start_time = time.time()

    file_count = split_file(input_file, lines_per_file)

    input_file_list = [f'split_inputs_reverse/inputrootlist_part{i}.txt' for i in range(file_count)]

## Method 1
    
#    run_macro(input_file_list, macro_file, output_prefix)

## Method 2
    with ThreadPoolExecutor(max_workers=len(input_file_list)) as executor:
        futures = []
        for i, input_file in enumerate(input_file_list):
            output_file = os.path.join(output_dir, f'selectedEvents_part{i}.root')
            futures.append(executor.submit(run_macro, input_file, output_file, macro_file))

        for future in futures:
            future.result()

'''
