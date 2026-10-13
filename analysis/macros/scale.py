import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import processor 
import hist
from coffea.util import load, save

def scale(filename):

	hists = load(filename)

	###
	# Rescaling MC histograms using the xsec weight
	###

	scale={}
	for dataset in hists['sumw'].keys():
		scale[dataset]=hists['sumw'][dataset]
		print('dataset: ', dataset)
	print('Sumw extracted')

	for key in hists.keys():
		if key=='sumw': continue
		for dataset in hists[key].keys():
			print('dataset in loop', dataset)
			if 'Mphi' in dataset:
				print(scale[dataset])
			#if 'Mphi' in dataset or  'Muon' in dataset or 'MET' in dataset or 'SingleElectron' in dataset or 'SinglePhoton' in dataset or 'EGamma' in dataset or 'BTagMu' in dataset: continue
			if  'Muon' in dataset or 'MET' in dataset or 'SingleElectron' in dataset or 'SinglePhoton' in dataset or 'EGamma' in dataset or 'BTagMu' in dataset: continue
			hists[key][dataset] *= 1/scale[dataset]
            #print('pass')
	print('Histograms scaled')


	###
	# Defining 'process', to aggregate different samples into a single process
	##
	sig_map = {}
	bkg_map = {}
	data_map = {}
	#bkg_map["VV"] = (["WW*","WZ*","ZZ*"],)
	bkg_map["WW"] = ["WW"]
	bkg_map["WZ"] = ["WZ"]
	bkg_map["ZZ"] = ["ZZ"]
	bkg_map["QCD Multijet"] = ["QCD"]
	#bkg_map["VV"] = ["_TuneCP5_13TeV-pythia8"]
	bkg_map["TT"] = ["TT"]
	#bkg_map["Single Top"] = ["ST"]
	bkg_map["ST"] = ["ST"]
	#bkg_map[r"W ($\ell\nu$) + Jets"]  = ["WJets"]
	#bkg_map[r"Z ($\ell\ell$) + Jets"] = ["DYJets"]
	bkg_map[r"Z ($\nu\nu$) + Jets"]   = ["Zto2Nu"]
	bkg_map[r"Z ($\ell\ell$) + Jets"] = ["DYto2L"]
	bkg_map[r"W ($\ell\nu$) + Jets"]  = ["WtoLNu"]
	#bkg_map[r"W ($\ell\nu$) + Jets PTLNu-200to400_1J"]  = ["WtoLNu-2Jets_PTLNu-200to400_1J"]
	#bkg_map[r"W ($\ell\nu$) 0J"]  = ["WtoLNu-2Jets_0J"]
	#bkg_map[r"W ($\ell\nu$) 1J"]  = ["WtoLNu-2Jets_1J"]
	#bkg_map[r"W ($\ell\nu$) 2J"]  = ["WtoLNu-2Jets_2J"]
	#bkg_map[r"W ($\ell\nu$) PTLNu"]  = ["WtoLNu-2Jets_PT"]
	bkg_map[r"G + Jets"]  = ["GJ_PTG"]
	#bkg_map[r"G + Jets PT-100to200"]  = ["GJ_PTG-100to200"]
	#bkg_map[r"G + Jets PT-200to400"]  = ["GJ_PTG-200to400"]
	#bkg_map[r"G + Jets PT-400to600"]  = ["GJ_PTG-400to600"]
	#bkg_map[r"G + Jets PT-600"]  = ["GJ_PTG-600"]
	#data_map["Muon"] = ["Muon"]
	data_map["MET"] = ["JetMET"]
	data_map["Muon"] = ["Muon"]
	#data_map["SingleElectron"] = ["SingleElectron"]
	#data_map["SinglePhoton"] = ["SinglePhoton"]
	data_map["EGamma"] = ["EGamma"]
	#data_map["BTagMu"] = ["BTagMu"]

	for signal in hists['sumw'].keys():
		if 'Mphi' not in signal: continue
		print(signal)
		#sig_map[signal] = signal  ## signals
		sig_map[signal] = [signal]  ## signals
	print('Processes defined')
	
	###
	# Storing signal and background histograms
	###
	bkg_hists={}
	sig_hists={}
	data_hists={}
	for key in hists.keys():
		bkg_hists[key]={}
		sig_hists[key]={}
		data_hists[key]={}
		for process in bkg_map.keys():
			for dataset in hists[key].keys():
				if not any(d in dataset for d in bkg_map[process]): continue
				#print('Adding',dataset,'to',process,'for variable',key)
				try:
					bkg_hists[key][process]+=hists[key][dataset]
				except:
					bkg_hists[key][process]=hists[key][dataset]
		for process in data_map.keys():
			for dataset in hists[key].keys():
				if not any(d in dataset for d in data_map[process]): continue
				#print('Adding',dataset,'to',process,'for variable',key)
				try:
					data_hists[key][process]+=hists[key][dataset]
				except:
					data_hists[key][process]=hists[key][dataset]
		for process in sig_map.keys():
			#print('Scaling '+ process +' by xsec '+str(xsec[process][2]))
			#sig_hists[key][process] *= xsec[str(process)][2]
			for dataset in hists[key].keys():
				#if not any(d in dataset for d in sig_map[process]): continue
				if dataset != process: continue
				#print('Adding',dataset,'to',process,'for variable',key)
				try:
					sig_hists[key][process]+=hists[key][dataset]
				except:
					sig_hists[key][process]=hists[key][dataset]
			#for signal in sig_hists[key]:
		
	print('Histograms grouped')

	return bkg_hists, sig_hists, data_hists

if __name__ == '__main__':
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('-f', '--file', help='file', dest='file')
	(options, args) = parser.parse_args()

	bkg_hists, sig_hists, data_hists = scale(options.file)
	name = options.file

	hists={
		'bkg': bkg_hists,
		'sig': sig_hists,
		'data': data_hists
	}
	save(hists,name.replace('.merged','.scaled'))
