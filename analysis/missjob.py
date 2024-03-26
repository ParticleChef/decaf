import fnmatch
import numpy as np
import numexpr
import subprocess
import os
import difflib
import glob
import sys
import json
import gzip
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-p', '--process', help='process', dest='process')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
parser.add_option('-d', '--dataset', help='see dataset you want', dest='dataset', default=False)
(options, args) = parser.parse_args()

metadata = options.metadata

with gzip.open('metadata/'+metadata+'.json.gz') as f:
    samplefiles = json.load(f)

totaljob = list(samplefiles.keys())
if options.dataset:
    f = open("./missJobs_"+str(options.process)+"_"+str(options.dataset)+".txt","w")
else:
    f = open("./missJobs_"+str(options.process)+".txt","w")

for job in samplefiles.keys():
    if options.dataset:
        if not str(options.dataset) in job: totaljob.remove(job)
#    f.write(job+"\n")

    dlist = glob.glob('./hists/'+options.process+'/'+job+'*')
    for i in dlist:
        if options.dataset:
            if not str(options.dataset) in i: continue
        ijob = i.split('/')[3].split('.')[0]
        totaljob.remove(ijob)
        #f.write(ijob+"\n")

for miss in totaljob:
    f.write(miss+"\n")
        
#cnt = 0        
#for d in lst:
#    totaljob = list()
#    for nj in range(njob[cnt]):
#        jobname = d + '____' + str(nj+1) + '_'
#        totaljob.append(jobname)
#
#    print(d, len(totaljob))
#    for miss in totaljob:
#        f.write(miss+"\n")
#    cnt = cnt+1
f.close()
