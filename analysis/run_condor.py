#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import numexpr
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset', default='')
parser.add_option('-e', '--exclude', help='exclude', dest='exclude', default='')
parser.add_option('-p', '--processor', help='processor', dest='processor', default='')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

os.system("mkdir -p hists/"+options.processor)

if options.tar:
	os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw_11_3_4.tgz '
			'--exclude=\'src/decaf/analysis/logs\' '
			'--exclude=\'src/decaf/analysis/plots\' '
			'--exclude=\'src/decaf/analysis/datacards\' '
			'--exclude=\'src/decaf/analysis/results\' '
			'--exclude=\'src/decaf/analysis/data/models\' '
			'--exclude=\'src/decaf/analysis/hists/*/*.futures\' '
			'--exclude=\'src/decaf/analysis/hists/*/*.merged\' '
			'--exclude=\'src/decaf/analysis/hists/*/*.reduced\' '
			'../../../../CMSSW_11_3_4')
	#os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../pylocal_3_8.tgz -C ~/.local/lib/python3.8/ site-packages')

if options.cluster == 'kisti':
	if options.copy:
		os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
		print('cmssw removed')
		os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
		os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/pylocal_3_8.tgz') 
		print('pylocal removed')
		os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal_3_8.tgz')
	jdl = """universe = vanilla
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh, /tmp/x509up_u556950957
Output = logs/condor/run/out/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/run/err/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/run/log/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(PROCESSOR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(PROCESSOR)/$ENV(SAMPLE).futures"
Arguments = $ENV(METADATA) $ENV(SAMPLE) $ENV(PROCESSOR) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(BTCN)
request_cpus = 8
request_memory = 7000
Queue 1"""

if options.cluster == 'lpc':
	if options.copy:
		os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
		os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal_3_8.tgz')
	jdl = """universe = vanilla
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh
Output = logs/condor/run/out/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/run/err/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/run/log/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(PROCESSOR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(PROCESSOR)/$ENV(SAMPLE).futures"
Arguments = $ENV(METADATA) $ENV(SAMPLE) $ENV(PROCESSOR) $ENV(CLUSTER) $ENV(USER) 
JobBatchName = $ENV(BTCN)
request_cpus = 8
request_memory = 16000
Queue 1"""

if options.cluster == 'knut3':
	if options.copy: # /pnfs/knu.ac.kr/data/cms/store/user/jhong/
		#os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://cluster142.knu.ac.kr//store/user/'+os.environ['USER']+'/cmssw_11_3_4.tgz')
		#os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://cluster142.knu.ac.kr//store/user/'+os.environ['USER']+'/pylocal_3_8.tgz')
		os.system('xrdcp -f ../../../../cmssw_11_3_4.tgz root://cluster142.knu.ac.kr//store/user/jhong/cmssw_11_3_4.tgz')
		os.system('xrdcp -f ../../../../pylocal_3_8.tgz root://cluster142.knu.ac.kr//store/user/jhong/pylocal_3_8.tgz')
	jdl = """universe = vanilla
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
use_x509userproxy = True
Transfer_Input_Files = run.sh
Output = logs/condor/run/out/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/run/err/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/run/log/$ENV(PROCESSOR)_$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(PROCESSOR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(PROCESSOR)/$ENV(SAMPLE).futures"
Arguments = $ENV(METADATA) $ENV(SAMPLE) $ENV(PROCESSOR) $ENV(CLUSTER) $ENV(USER) 
run_as_owner = True
JobBatchName = $ENV(BTCN)
request_memory = 2048 MB
Queue 1"""

jdl_file = open("run.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

with gzip.open('metadata/'+options.metadata+'.json.gz') as fin:
	datadef = json.load(fin)

for dataset, info in datadef.items():
	if options.dataset:
		if not any(_dataset in dataset for _dataset in options.dataset.split(',')): continue
	if options.exclude:
		if any(_dataset in dataset for _dataset in options.exclude.split(',')): continue
	os.system('mkdir -p logs/condor/run/err/')
	os.system('rm -rf logs/condor/run/err/*'+options.processor+'*'+dataset+'*')
	os.system('mkdir -p logs/condor/run/log/')
	os.system('rm -rf logs/condor/run/log/*'+options.processor+'*'+dataset+'*')
	os.system('mkdir -p logs/condor/run/out/')
	os.system('rm -rf logs/condor/run/out/*'+options.processor+'*'+dataset+'*')
	os.environ['SAMPLE'] = dataset
	os.environ['BTCN'] = dataset.split('____')[0]
	os.environ['PROCESSOR'] = options.processor
	os.environ['METADATA'] = options.metadata
	os.environ['CLUSTER'] = options.cluster
	os.system('condor_submit run.submit')
os.system('rm run.submit')
