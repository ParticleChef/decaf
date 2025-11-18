#!/usr/bin/env python
import logging
import numpy as np
import awkward as ak
import json
import copy
from collections import defaultdict
from coffea import processor
import cachetools
import hist
from coffea.analysis_tools import Weights, PackedSelection
from coffea.lumi_tools import LumiMask
from coffea.util import load, save
from optparse import OptionParser
from coffea.nanoevents.methods import vector
import gzip

def update(events, collections):
	"""Return a shallow copy of events array with some collections swapped out"""
	out = events
	for name, value in collections.items():
		out = ak.with_field(out, value, name)
	return out

class AnalysisProcessor(processor.ProcessorABC):

	lumis = { 
		#https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiRecommendationsRun3
		'2022pre': 7.980,
		'2022post': 26.67,
	}

	lumiMasks = {
		'2022pre': LumiMask("data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json"),
	}
	
	met_filters = {
	# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC
		'2022pre': [
				'goodVertices', 
				'globalSuperTightHalo2016Filter', 
				'EcalDeadCellTriggerPrimitiveFilter', 
				'BadPFMuonFilter', 
				'BadPFMuonDzFilter',
				'hfNoisyHitsFilter',
				'eeBadScFilter', 
				]
	}
			
	def __init__(self, year, xsec, corrections, ids, common):

		self._year = year
		self._lumi = 1000.*float(AnalysisProcessor.lumis[year])
		self._xsec = xsec
		self._systematics = False
		self._skipJER = True

		self._samples = {
			'total' :('GJ','JetMET'),
			'pho200' :('GJ', 'JetMET'),
		}
		

		self._met_triggers = { ## NOT UPDATED for Run3
			'2022pre': [
				'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
				'PFMETNoMu120_PFMHTNoMu120_IDTight'
			]
		}
		self._photon_ref_triggers = {
			'2022pre':
				[
				'PFHT1050'
			]
		}
		self._photon_triggers = {
			'2022pre':
				[
				'Photon200'
			]
		}
		self._singlemuon_triggers = {
			'2022pre':
				[
				'Mu50',
				'CascadeMu100',
				'HighPtTkMu100'
			]
		}
		self._corrections = corrections
		self._ids = ids
		self._common = common

		self.make_output = lambda: {
			'sumw': 0.,
			'template': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(100,0,1000, name='phopt', label='Photon Pt'),
                #hist.axis.Variable([0,1], name='isPass', label='Photon200'),
				#hist.axis.Variable([40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300], name='fjmass', label=r'AK15 Jet $m_{sd}$'),
				storage=hist.storage.Weight(),
			),
			'phopt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='phopt', label='Leading Photon Pt'),
				storage=hist.storage.Weight(),
			),
			'phoeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='phoeta', label='Leading Photon Eta'),
				storage=hist.storage.Weight(),
			),
			'phophi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='phophi', label='Leading Photon Phi'),
				storage=hist.storage.Weight(),
			),
	}

	def process(self, events):
		isData = not hasattr(events, "genWeight")
		if isData:
			# Nominal JEC are already applied in data
			return self.process_shift(events, None)
		
		if not isData:
			return self.process_shift(events, None)
	
	def process_shift(self, events, shift_name):

		dataset = events.metadata['dataset']

		selected_regions = []
		for region, samples in self._samples.items():
			for sample in samples:
				if sample not in dataset: continue
				selected_regions.append(region)
		
		isData = not hasattr(events, "genWeight")
		
		################
		## Initialize selection, weights and output histogram
		################
		selection = PackedSelection(dtype="uint64")
		weights = Weights(len(events), storeIndividual=True)
		output = self.make_output()

		
		if shift_name is None and not isData:
			output['sumw'] = ak.sum(events.genWeight)

		###
		#Getting corrections, ids from .coffea files
		###

		isLooseMuon	 = self._ids['isLooseMuon']	 
		isTightMuon	 = self._ids['isTightMuon']	 
		#isLooseElectron = self._ids['isLooseElectron'] 
		#isTightElectron = self._ids['isTightElectron'] 
		isLoosePhoton   = self._ids['isLoosePhoton']   
		isTightPhoton   = self._ids['isTightPhoton']   

		###
		#Initialize global quantities (MET ecc.)
		###

		npv = events.PV.npvsGood
		met = events.MET

		###
		#Initialize physics objects
		###

		pho = events.Photon
		pho['isloose'] = isLoosePhoton(pho,self._year)
		pho['istight'] = isTightPhoton(pho,self._year)

		pho_loose=pho[pho.isloose]
		pho_tight=pho[pho.istight]
		pho_ntot = ak.num(pho, axis=1)
		pho_nloose = ak.num(pho_loose, axis=1)
		pho_ntight = ak.num(pho_tight, axis=1)
		
		# define leading pho
		leading_pho = ak.firsts(pho_tight)
		
		###
		# Selections
		###

		lumimask = np.ones(len(events), dtype='bool')
		if isData:
			lumimask = AnalysisProcessor.lumiMasks[self._year](events.run, events.luminosityBlock)
		selection.add('lumimask', lumimask)

		met_filters =  np.ones(len(events), dtype='bool')
		for flag in AnalysisProcessor.met_filters[self._year]:
			met_filters = met_filters & events.Flag[flag]
		selection.add('met_filters',met_filters)

		triggers = np.zeros(len(events), dtype='bool')
		for path in self._met_triggers[self._year]:
			if not hasattr(events.HLT, path): continue
			triggers = triggers | events.HLT[path]
		#selection.add('met_triggers', triggers)

        ## Reference trigger
		triggers = np.zeros(len(events), dtype='bool')
		for path in self._photon_ref_triggers[self._year]:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('reference_triggers', ak.to_numpy(triggers))

        ## Photon Trigger
		triggers = np.zeros(len(events), dtype='bool')
		for path in self._photon_triggers[self._year]:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('single_photon_triggers', ak.to_numpy(triggers))

		

		selection.add('photon_selection', ((pho_ntight==1)&(pho_nloose==1)))


		regions = {
			'total': [
					'lumimask',
					'met_filters',
					'photon_selection',
                    'reference_triggers'
			],
			'pho200': [
					'lumimask',
					'met_filters',
					'photon_selection',
                    'reference_triggers',
					'single_photon_triggers'
			],
		}

		def normalize(val, cut):
			if cut is None:
				ar = ak.to_numpy(ak.fill_none(val, np.nan))
				return ar
			else:
				ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
				return ar
				
		def fill(region, systematic):
			cut = selection.all(*regions[region])
			#sname = 'nominal' if systematic is None else systematic
			#if systematic in weights.variations:
			#	weight = weights.weight(modifier=systematic)[cut]
			#else:
			#	weight = weights.weight()[cut]
			output['template'].fill(
				  region=region,
				  #systematic=sname,
                  phopt=normalize(leading_pho.pt, cut),
#				  weight=weight
			)
			if systematic is None:
				variables = {
					'phopt':				  leading_pho.pt,
					'phophi':				 leading_pho.phi,
					'phoeta':				 leading_pho.eta,
				}
				for variable in output:
					if variable not in variables:
						continue
					normalized_variable = {variable: normalize(variables[variable],cut)}
					output[variable].fill(
						region=region,
						**normalized_variable,
			#			weight=weight,
					)


		if shift_name is None:
			systematics = [None] + list(weights.variations)
		else:
			systematics = [shift_name]
			
		for region in regions:
			if region not in selected_regions: continue


			###
			# Adding recoil and minDPhi requirements
			###

			#if True:
			#	selection.add('test_selection', np.ones_like(leading_mu.pt))
			#	
			#	if region == 'nocr':
			#		regions[region].insert(2, 'test_selection')
			#	else:
			#		regions[region].insert(3, 'test_selection')

			for systematic in systematics:
				if isData and systematic is not None:
					continue
				fill(region, systematic)


		scale = 1
		if self._xsec[dataset]!= -1: 
			scale = self._lumi*self._xsec[dataset]

		for key in output:
			if key=='sumw': 
				continue
			output[key] *= scale #* skimrate[dataset.split('____')]
				
		return output

	def postprocess(self, accumulator):

		return accumulator

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-y', '--year', help='year', dest='year')
	parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
	parser.add_option('-n', '--name', help='name', dest='name')
	(options, args) = parser.parse_args()


	with gzip.open('metadata/'+options.metadata+'.json.gz') as fin:
		samplefiles = json.load(fin)
		xsec = {k: v['xs'] for k,v in samplefiles.items()}

	corrections = load('data/corrections.coffea')
	ids		 = load('data/ids.coffea')
	common	  = load('data/common.coffea')

	processor_instance=AnalysisProcessor(year=options.year,
										 xsec=xsec,
										 corrections=corrections,
										 ids=ids,
										 common=common)

	save(processor_instance, 'data/hadmonotop'+options.name+'.processor')
