import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import processor 
import hist
from coffea.util import load, save

xsec = {
	'sig_Mphi-1000_Mchi-10': ('sig_Mphi-1000_Mchi-10', 'MC', 0.4298),
	'sig_Mphi-1000_Mchi-1000': ('sig_Mphi-1000_Mchi-1000', 'MC', 6.641e-05),
	'sig_Mphi-1000_Mchi-150': ('sig_Mphi-1000_Mchi-150', 'MC', 0.4231),
	'sig_Mphi-1000_Mchi-50': ('sig_Mphi-1000_Mchi-50', 'MC', 0.4283),
	'sig_Mphi-100_Mchi-150': ('sig_Mphi-100_Mchi-150', 'MC', 0.1315),
	'sig_Mphi-1245_Mchi-625': ('sig_Mphi-1245_Mchi-625', 'MC', 0.03257),
	'sig_Mphi-1250_Mchi-150': ('sig_Mphi-1250_Mchi-150', 'MC', 0.1663),
	'sig_Mphi-1495_Mchi-750': ('sig_Mphi-1495_Mchi-750', 'MC', 0.01422),
	'sig_Mphi-1500_Mchi-150': ('sig_Mphi-1500_Mchi-150', 'MC', 0.07259),
	'sig_Mphi-150_Mchi-150': ('sig_Mphi-150_Mchi-150', 'MC', 0.2064),
	'sig_Mphi-1700_Mchi-800': ('sig_Mphi-1700_Mchi-800', 'MC', 0.02241),
	'sig_Mphi-1750_Mchi-150': ('sig_Mphi-1750_Mchi-150', 'MC', 0.03415),
	'sig_Mphi-1750_Mchi-700': ('sig_Mphi-1750_Mchi-700', 'MC', 0.0269),
	'sig_Mphi-195_Mchi-100': ('sig_Mphi-195_Mchi-100', 'MC', 5.346),
	'sig_Mphi-1995_Mchi-1000': ('sig_Mphi-1995_Mchi-1000', 'MC', 0.003229),
	'sig_Mphi-2000_Mchi-150': ('sig_Mphi-2000_Mchi-150', 'MC', 0.01696),
	'sig_Mphi-2000_Mchi-1500': ('sig_Mphi-2000_Mchi-1500', 'MC', 5.795e-06),
	'sig_Mphi-2000_Mchi-500': ('sig_Mphi-2000_Mchi-500', 'MC', 0.01562),
	'sig_Mphi-200_Mchi-10': ('sig_Mphi-200_Mchi-10', 'MC', 63.67),
	'sig_Mphi-200_Mchi-150': ('sig_Mphi-200_Mchi-150', 'MC', 0.2782),
	'sig_Mphi-200_Mchi-20': ('sig_Mphi-200_Mchi-20', 'MC', 64.21),
	'sig_Mphi-200_Mchi-30': ('sig_Mphi-200_Mchi-30', 'MC', 63.39),
	'sig_Mphi-200_Mchi-40': ('sig_Mphi-200_Mchi-40', 'MC', 63.79),
	'sig_Mphi-200_Mchi-50': ('sig_Mphi-200_Mchi-50', 'MC', 63.16),
	'sig_Mphi-2245_Mchi-1125': ('sig_Mphi-2245_Mchi-1125', 'MC', 0.001634),
	'sig_Mphi-2250_Mchi-150': ('sig_Mphi-2250_Mchi-150', 'MC', 0.008869),
	'sig_Mphi-2495_Mchi-1250': ('sig_Mphi-2495_Mchi-1250', 'MC', 0.0008546),
	'sig_Mphi-2500_Mchi-150': ('sig_Mphi-2500_Mchi-150', 'MC', 0.004808),
	'sig_Mphi-2500_Mchi-2000': ('sig_Mphi-2500_Mchi-2000', 'MC', 3.947e-07),
	'sig_Mphi-2500_Mchi-750': ('sig_Mphi-2500_Mchi-750', 'MC', 0.004033),
	'sig_Mphi-295_Mchi-150': ('sig_Mphi-295_Mchi-150', 'MC', 2.595),
	'sig_Mphi-2995_Mchi-1500': ('sig_Mphi-2995_Mchi-1500', 'MC', 0.0002481),
	'sig_Mphi-3000_Mchi-1000': ('sig_Mphi-3000_Mchi-1000', 'MC', 0.00115),
	'sig_Mphi-3000_Mchi-2000': ('sig_Mphi-3000_Mchi-2000', 'MC', 6.318e-07),
	'sig_Mphi-300_Mchi-100': ('sig_Mphi-300_Mchi-100', 'MC', 20.63),
	'sig_Mphi-300_Mchi-300': ('sig_Mphi-300_Mchi-300', 'MC', 0.02127),
	'sig_Mphi-400_Mchi-150': ('sig_Mphi-400_Mchi-150', 'MC', 8.83),
	'sig_Mphi-495_Mchi-250': ('sig_Mphi-495_Mchi-250', 'MC', 0.7786),
	'sig_Mphi-500_Mchi-150': ('sig_Mphi-500_Mchi-150', 'MC', 4.759),
	'sig_Mphi-500_Mchi-500': ('sig_Mphi-500_Mchi-500', 'MC', 0.002621),
	'sig_Mphi-625_Mchi-150': ('sig_Mphi-625_Mchi-150', 'MC', 2.356),
	'sig_Mphi-745_Mchi-325': ('sig_Mphi-745_Mchi-325', 'MC', 1.009),
	'sig_Mphi-750_Mchi-150': ('sig_Mphi-750_Mchi-150', 'MC', 1.26),
	'sig_Mphi-995_Mchi-500': ('sig_Mphi-995_Mchi-500', 'MC', 0.08102),
}
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
