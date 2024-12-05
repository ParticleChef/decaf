import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import time

def split_file(input_file, lines_per_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    total_lines = len(lines)

    ## How many jobs
    file_count = (total_lines // lines_per_file) + (1 if total_lines % lines_per_file !=0 else 0)

    os.makedirs('split_inputs', exist_ok=True)
    for i in range(file_count):
        #with open(f'split_inputs/inputrootlist_part_{i}.txt', 'w') as f_part:
        with open(f'split_inputs/evt_ids_egm_ntuple_part_{i}.txt', 'w') as f_part:
            for line in lines[i*lines_per_file:(i+1)*lines_per_file]:
                f_part.write(line)
    return file_count

def run_macro(input_file, output_file, macro_file):
#    for i, input_file in enumerate(input_file_list):
    #cmd = ['root', '-l', '-b', '-q', f'{macro_file}("{input_file}", "{output_file}", "Events", "EventNo_mTop.txt")']
    cmd = ['root', '-l', '-b', '-q', f'{macro_file}("{input_file}", "{output_file}", "Events", "event_ids_2018_egm_from_ntuple.txt")']
    subprocess.run(cmd)

if __name__ == "__main__":
    input_file =  'onefileegm.txt'  ## all root files 
    #input_file = 'inputrootlist.txt'   ## part of root files
    #input_file = 'inputrootMichael.txt'   ## Michael's ntuple
    lines_per_file = 35
    macro_file = 'test3diff.C'    #### test2.C: compare event number, test3.C: compare three components

    #output_dir = 'result_files_again'
    #output_dir = '/10T/scratch/jhong/selectedEvents_result_files'
    output_dir = './output' #'/10T/scratch/jhong/onlyMoritz'
    os.makedirs(output_dir, exist_ok=True)

    start_time = time.time()

    file_count = split_file(input_file, lines_per_file)

    input_file_list = [f'split_inputs/evt_ids_egm_ntuple_part_{i}.txt' for i in range(file_count)]

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

print('Total time: ', str(time.time()-start_time), ' seconds')
