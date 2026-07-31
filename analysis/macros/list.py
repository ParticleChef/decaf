#!/usr/bin/env python
import os
import sys
import uproot
from data.process import *
from optparse import OptionParser
import json
import gzip

parser = OptionParser()
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
parser.add_option('-p', '--pack', help='pack', dest='pack')
parser.add_option('-s', '--special', help='special', dest='special')
parser.add_option('-c', '--custom', action='store_true', dest='custom')
parser.add_option('-k', '--skip', help='skip', dest='skip')
parser.add_option('-r', '--remove', action='store_true', dest='remove')
(options, args) = parser.parse_args()

globalredirect = "root://xrootd-cms.infn.it/"
campaigns ={}
campaigns['2016preVFP'] = '*UL*16preVFP*JMENano'
campaigns['2016postVFP'] = '*UL*16postVFP*JMENano'
campaigns['2017'] = '*UL*17*JMENano'
campaigns['2018'] = '*UL*18*JMENano'
campaigns['2022pre'] = '*Run3*22Nano*v2'
campaigns['2023pre'] = '*Run3*22Nano*v2'

eos = "root://dcache-cms-xrootd.desy.de:1094/"
custom={}
custom['2016preVFP'] = ["/store/user/nshadski/customNano",
				"/store/user/empfeffe/customNano",
				"/store/user/momolch/customNano",
				"/store/user/swieland/customNano",
				"/store/user/mwassmer/customNano"]

custom['2016postVFP'] = ["/store/user/nshadski/customNano",
				"/store/user/empfeffe/customNano",
				"/store/user/momolch/customNano",
				"/store/user/swieland/customNano",
				"/store/user/mwassmer/customNano"]

custom['2017'] = ["/store/user/swieland/customNano",
				"/store/user/momolch/customNano",
				"/store/user/mwassmer/customNano"]

custom['2018'] = ["/store/user/mwassmer/customNano",
				"/store/user/swieland/customNano"]
custom['2023pre'] = ["/store/group/lpcmetx/Monotop/NanoAOD/2023",
				"/store/group/lpcmetx/Monotop/NanoAOD/2023"]
custom['2023post'] = ["/store/user/jhong/monotopRun3/2023BPix/NanoAOD"]


def split(arr, size):
	arrs = []
	while len(arr) > size:
		pice = arr[:size]
		arrs.append(pice)
		arr = arr[size:]
	arrs.append(arr)
	return arrs

def find(_list):
	if not _list:
		return []
	files=[]
	print('Looking into',_list)
	for path in _list:
		command='xrdfs '+eos+' ls '+path
		results=os.popen(command).read()
		files.extend(results.split())
	if not any('.root' in _file for _file in files):
		files=find(files)
	return files

xsections={}
for k,v in processes.items():
	if v[1]=='MC':
		if not isinstance(k, str):
			if options.year!=str(k[1]): continue
			xsections[k[0]] = v[2]
		else: 
			xsections[k] = v[2]
	else:
		xsections[k] = -1

if options.skip:
	try:
		os.system('ls '+options.skip)
	except:
		sys.exit('File',options.skip,'does not exist')
		
	skip = []
	corrupted = open(options.skip, 'r')
	for rootfile in corrupted.readlines():
		skip.append(rootfile.strip().split('store')[1])

removed = []
datadef = {}
datasets = []
for dataset in xsections.keys():
	print("dataset: ", dataset)
	if options.dataset:
		if not any(_dataset in dataset for _dataset in options.dataset.split(',')): continue
	xs = xsections[dataset]
	if options.custom:
		redirect = globalredirect#eos
		urllist = []
		for folder in custom[options.year]:
			path=folder+'/'+dataset
			urllist += find([path])
		for url in urllist[:].copy():
			print(url)
			if options.year not in url:
					urllist.remove(url)
					continue
			if 'Data' in url and 'KITv2' in url:
					urllist.remove(url)
					continue
			if 'failed' in url: 
					urllist.remove(url)
					continue
			if '.root' not in url: 
					urllist.remove(url)
					continue
			if options.skip and url.split('store')[-1] in skip:
					urllist.remove(url)
					print(url,'found in',options.skip)
					continue
			if options.remove:
					try:
						infile = uproot.open(redirect+url)
					except:
						print("File",redirect+url,"is corrupted, removing.")
						urllist.remove(url)
						removed.append(url)
						continue
					else:
						del infile

	else:
		redirect = globalredirect
		print("Searching for",dataset,"in centrally produced NanoAOD")
		if 'Muon' in dataset:
			query="dasgoclient --query=\"dataset dataset=/"+dataset+"/Run2022C*/NANOAOD*\""
		else:
			query="dasgoclient --query=\"dataset dataset=/"+dataset+"/"+campaigns[options.year]+"*/NANOAOD*\""
		dataset=os.popen(query).read().split("\n")[0]
		print('Dataset is:', dataset)
		query="dasgoclient --query=\"file dataset="+dataset+"\""
		urllist = os.popen(query).read().split("\n")
	for url in urllist[:]:
		urllist[urllist.index(url)]=redirect+url
	print('list lenght:',len(urllist))
	if options.special:
		for special in options.special.split(','):
			sdataset, spack = special.split(':')
			if sdataset in dataset:
				print('Packing',spack,'files for dataset',dataset)
				urllists = split(urllist, int(spack))
			else:
				print('Packing',int(options.pack),'files for dataset',dataset)
				urllists = split(urllist, int(options.pack))
	else:
		print('Packing',int(options.pack),'files for dataset',dataset)
		urllists = split(urllist, int(options.pack))
	print(len(urllists))
	if urllist:
		datasetname = dataset.split("/")[1]
		for i in range(0,len(urllists)) :
			datadef[datasetname+"____"+str(i+1)+"_"] = {
				'files': urllists[i],
				'xs': xs,
			}
		
json_output = "metadata/"+options.metadata+".json.gz"
with gzip.open(json_output, "wt") as fout:
	json.dump(datadef, fout, indent=4)
	#fp.write("\n")

if options.remove:
	lists = "data/removed_files.txt"
	with open(lists, "w") as fout:
		fout.writelines(removed)
