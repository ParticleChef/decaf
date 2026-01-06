import os
import json
import time
from optparse import OptionParser
import gzip
# 
parser = OptionParser()
parser.add_option('-p', '--processor', help='processor', dest='processor')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-w', '--workers', help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=1)
(options, args) = parser.parse_args()

start_time = time.time()

with gzip.open("metadata/"+options.metadata+".json.gz") as fin:
    samplefiles = json.load(fin)
i, idx = 0, 0

for dataset, info in samplefiles.items():
    if options.dataset and options.dataset not in dataset: continue
    # check players in lane
    #print('total samples: '+str(total_samples))
    sleep_time = 10
    stream = os.popen(r"ps -eo ppid=,args= | awk '$1==1 && $0 ~ /python3 run\.py/ {c++} END{print c+0}'")
    if int(stream.read()) < 50:
        # nohup job
        #print('ds')
        os.system('nohup python3 run.py -p '+options.processor+' -m '+options.metadata+' -d '+dataset+' -w '+str(options.workers)+' > log/'+dataset+'.log &')
        idx += 1
        i = 0
        time.sleep(1)
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
            time.sleep(sleep_time)
            stream = os.popen(r"ps -eo ppid=,args= | awk '$1==1 && $0 ~ /python3 run\.py/ {c++} END{print c+0}'")
            if int(stream.read()) < 50:
                break
        os.system('nohup python3 run.py -p '+options.processor+' -m '+options.metadata+' -d '+dataset+' -w '+str(options.workers)+' > log/'+dataset+'.log &')
        idx += 1
        i = 0
        time.sleep(1)
        continue
        #print("this!")

print('total '+str(idx)+' jobs are submitted!')
print('all jobs are submitted!')
print('starting monitoring...')
while True:
    if options.dataset: 
        stream = os.popen("ls -lh hists/"+options.processor+"| grep "+options.dataset+" | wc -l")
    else:
        stream = os.popen("ls -lh hists/"+options.processor+"| grep .futures |wc -l")
    print('----------------------------------------------------------')
    print('For now "'+str(stream.read())+'"       players has finished')
    stream = os.popen(r"ps -eo ppid=,args= | awk '$1==1 && $0 ~ /python3 run\.py/ {c++} END{print c+0}'")
    print('Still '+str(int(stream.read())/2)+' jobs are running')
    print('sleeping for '+str(sleep_time*i)+' seconds...')
    print('----------------------------------------------------------')
    time.sleep(sleep_time)
    i += 1
    stream = os.popen("ps -ef | grep "+options.processor+" | wc -l")
    if int(stream.read()) < 4:
        break

if options.dataset:
    os.system('python3 missjob.py -p '+options.processor+' -m '+options.metadata+' -d '+options.dataset)
    missjobs = open("./missJobs_"+str(options.processor)+"_"+str(options.dataset)+".txt", "r")
else:
    os.system('python3 missjob.py -p '+options.processor+' -m '+options.metadata)
    missjobs = open("./missJobs_"+str(options.processor)+".txt", "r")

missjob = missjobs.readlines()

if len(missjob) < 2:
    print('No miss job! Yeah!')
else:
    print('Oh... Check the missJobs_'+str(options.processor)+".txt")
#
#while True:
#    if options.dataset:
#        os.system('python3 missjob.py -p '+options.processor+' -m '+options.metadata+' -d '+options.dataset)
#        missjobs = open("./missJobs_"+str(options.processor)+"_"+str(options.dataset)+".txt", "r")
#    else:
#        os.system('python3 missjob.py -p '+options.processor+' -m '+options.metadata)
#        missjobs = open("./missJobs_"+str(options.processor)+".txt", "r")
#
#    missjob = missjobs.readlines()
#
#    if len(missjob) == 0:
#        print('No miss job! Yeah!')
#        break
#    
#    
#    for dataset in missjob:
#        sleep_time = 10
#        stream = os.popen("ps -ef | grep 'python3 run.py' | wc -l")
#        if int(stream.read()) < 40:
#            # nohup job
#            os.system('nohup python3 run.py -p '+options.processor+' -m '+options.metadata+' -d '+dataset+' -w '+str(options.workers)+' > log/'+missjob+'.log &')
#            idx += 1
#            i = 0
#            time.sleep(1)
#            continue
#        else:
#            while True:
#                stream = os.popen("ps -ef | grep 'python3 run.py' | wc -l")
#                print('----------------------------------------------------------')
#                print('now "'+str(stream.read())+'"       players in lane')
#                print('total '+str(idx)+' jobs are submitted')
#                print('sleeping for '+str(sleep_time*i)+' seconds...')
#                print('----------------------------------------------------------')
#                time.sleep(sleep_time)
#                i += 1
#                time.sleep(sleep_time)
#                stream = os.popen("ps -ef | grep 'python3 run.py' | wc -l")
#                if int(stream.read()) < 40:
#                    break
#            os.system('nohup python3 run.py -p '+options.processor+' -m '+options.metadata+' -d '+dataset+' -w '+str(options.workers)+' > log/'+missjob+'.log &')
#            idx += 1
#            i = 0
#            time.sleep(1)
#            continue
#    print('starting monitoring...')
#
#    while True:
#        stream = os.popen("ps -ef | grep 'python3 run.py' | wc -l")
#        print('----------------------------------------------------------')
#        print('Still '+int(stream.read())+' jobs are running')
#        print('sleeping for '+str(sleep_time*i)+' seconds...')
#        print('----------------------------------------------------------')
#        time.sleep(sleep_time)
#        stream = os.popen("ps -ef | grep '"+options.processor+"' | wc -l")
#        if int(stream.read()) < 4:
#            break

print('Job done! Taiwoo is happy!')
# get metadata
print('Total time: ', str(time.time()-start_time), ' seconds')
