import os
import json
import time
from optparse import OptionParser
import gzip

start_time = time.time()
# 
parser = OptionParser()
parser.add_option('-p', '--processor', help='processor', dest='processor')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-w', '--workers', help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=1)
(options, args) = parser.parse_args()

with gzip.open("metadata/"+options.metadata+".json.gz") as fin:
    samplefiles = json.load(fin)
i, idx = 0, 0

for dataset, info in samplefiles.items():
    if options.dataset and options.dataset not in dataset: continue
    # check players in lane
    #print('total samples: '+str(total_samples))
    sleep_time = 10
    stream = os.popen(r"ps -eo ppid=,args= | awk '$1==1 && $0 ~ /python3 run\.py/ {c++} END{print c+0}'")
    if int(stream.read()) < 100:
        # nohup job
        #print('ds')
        print("sending", dataset)
        os.system('nohup python3 run.py -p '+options.processor+' -m '+options.metadata+' -d '+dataset+' -w '+str(options.workers)+' > log/'+dataset+'.log &')
        idx += 1
        i = 0
        time.sleep(0.5)
        continue
    else:
        while True:
            stream = os.popen(r"ps -eo ppid=,args= | awk '$1==1 && $0 ~ /python3 run\.py/ {c++} END{print c+0}'")
            print('----------------------------------------------------------')
            print('now "'+str(stream.read())+'"       players in lane')
            print('total '+str(idx)+' jobs are submitted')
            print('sleeping for '+str(sleep_time*i)+' seconds...')
            print('----------------------------------------------------------')
            time.sleep(sleep_time)
            i += 1
            stream = os.popen(r"ps -eo ppid=,args= | awk '$1==1 && $0 ~ /python3 run\.py/ {c++} END{print c+0}'")
            if int(stream.read()) < 100:
                break
        os.system('nohup python3 run.py -p '+options.processor+' -m '+options.metadata+' -d '+dataset+' -w '+str(options.workers)+' > log/'+dataset+'.log &')
        idx += 1
        i = 0
        continue

print('total '+str(idx)+' jobs are submitted!')
print('Total time: '+str(time.time()-start_time)+' seconds')
print('all jobs are submitted!')
print('starting monitoring...')
while True:
    stream = os.popen("ls -lh hists/"+options.processor+"| grep "+options.dataset+" | wc -l")
    print('----------------------------------------------------------')
    print('For now "'+str(stream.read())+'"       players has finished')
    stream = os.popen("ls -lh hists/"+options.processor+"| grep "+options.dataset+" | wc -l")
    print('Still '+str(idx-int(stream.read()))+' jobs are running')
    print('sleeping for '+str(sleep_time*i)+' seconds...')
    print('----------------------------------------------------------')
    time.sleep(sleep_time)
    i += 1
    stream = os.popen("ls -lh hists/"+options.processor+"| grep "+options.dataset+" | wc -l")
    if int(stream.read()) == idx:
        break
print('Job done! Taiwoo is happy!')
print('Total time: '+str(time.time()-start_time)+' seconds')