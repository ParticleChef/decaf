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
		'2022pre' : LumiMask("data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json"),
		'2022post' : LumiMask("data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json"),
		'2023pre': LumiMask("data/lumiMask/Cert_Collisions2023_366442_370790_Golden.json"),
		'2023post': LumiMask("data/lumiMask/Cert_Collisions2023_366442_370790_Golden.json"),
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
			'basic' :('GJ','JetMET'),
			'oneA' :('GJ','JetMET'),
			'eff' :('GJ', 'JetMET'),
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
		self._corrections = corrections
		self._ids = ids
		self._common = common

		self.make_output = lambda: {
			'sumw': 0.,
			'template': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='phopt', label='Photon Pt'),
				hist.axis.Variable([-0.5,0.5,1.5], name='PFHT1050', label='Reference trigger'),#growth=False),
				hist.axis.Variable([-0.5,0.5,1.5], name='Photon200',label='Photon trigger'),#growth=False),
				#hist.axis.Variable([40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300], name='fjmass', label=r'AK15 Jet $m_{sd}$'),
				storage=hist.storage.Weight(),
			),
			'phopt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(24,0,1200, name='phopt', label='Leading Photon Pt'),
				hist.axis.Integer(0,2, name='PFHT1050', label='Reference trigger'),
				hist.axis.Integer(0,2, name='Photon200', label='Photon trigger'),
				storage=hist.storage.Weight(),
			),
			'phoeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='phoeta', label='Leading Photon Eta'),
				hist.axis.Integer(0,2, name='PFHT1050', label='Reference trigger'),
				hist.axis.Integer(0,2, name='Photon200', label='Photon trigger'),
				storage=hist.storage.Weight(),
			),
			'phophi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='phophi', label='Leading Photon Phi'),
				hist.axis.Integer(0,2, name='PFHT1050', label='Reference trigger'),
				hist.axis.Integer(0,2, name='Photon200', label='Photon trigger'),
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
		isLooseElectron = self._ids['isLooseElectron'] 
		isTightElectron = self._ids['isTightElectron'] 
		isLoosePhoton   = self._ids['isLoosePhoton']   
		isTightPhoton   = self._ids['isTightPhoton']   
		isGoodAK4	   = self._ids['isGoodAK4']

		###
		#Initialize global quantities (MET ecc.)
		###

		npv = events.PV.npvsGood
		met = events.MET

		###
		#Initialize physics objects
		###

		mu = events.Muon
		mu['isloose'] = isLooseMuon(mu,self._year)
		mu_loose = mu[mu.isloose]

		e = events.Electron
		e['isloose'] = isLooseElectron(e,self._year)
		e_loose = e[e.isloose]

		pho = events.Photon
		pho['isloose'] = isLoosePhoton(pho,self._year)
		pho['istight'] = isTightPhoton(pho,self._year)

		pho_loose=pho[pho.isloose]
		pho_tight=pho[pho.istight]
		pho_ntot = ak.num(pho, axis=1)
		pho_nloose = ak.num(pho_loose, axis=1)
		pho_ntight = ak.num(pho_tight, axis=1)
		leading_pho = ak.firsts(pho_tight)
		
		j = events.Jet
		j['isgood'] = isGoodAK4(j, self._year)
		j['isclean'] = (
			ak.all(j.metric_table(mu_loose) > 0.4, axis=2)
			& ak.all(j.metric_table(e_loose) > 0.4, axis=2)
		)
		j_good = j[j.isgood]
		j_clean = j_good[j_good.isclean]
		j_ht = ak.sum(j_clean.pt, axis=1)

		
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

		

		selection.add('jet_ht1500', (j_ht > 1500))
		selection.add('photon_selection', ((pho_ntight==1)&(pho_nloose==1)))


		regions = {
			'basic': [
					'lumimask',
					'met_filters',
					#'jet_ht1500',
					#'photon_selection'
			],
			'oneA': [
					'lumimask',
					'met_filters',
					'photon_selection',
					#'jet_ht1500',
					#'reference_triggers'
			],
			'eff': [
					'lumimask',
					'met_filters',
					'photon_selection',
					'jet_ht1500',
					#'reference_triggers',
					#'single_photon_triggers'
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
			output['template'].fill(
				  region=region,
				  PFHT1050 =normalize(events.HLT.PFHT1050,cut),
				  Photon200=normalize(events.HLT.Photon200,cut),
				  phopt=normalize(leading_pho.pt, cut),
#				  weight=weight
			)
			output['phopt'].fill(
				  region=region,
				  phopt=normalize(leading_pho.pt, cut),
				  PFHT1050 =normalize(events.HLT.PFHT1050,cut),
				  Photon200=normalize(events.HLT.Photon200,cut),
#				  weight=weight
			)
			output['phophi'].fill(
				  region=region,
				  phophi=normalize(leading_pho.phi, cut),
				  PFHT1050 =normalize(events.HLT.PFHT1050,cut),
				  Photon200=normalize(events.HLT.Photon200,cut),
#				  weight=weight
			)
			output['phoeta'].fill(
				  region=region,
				  phoeta=normalize(leading_pho.eta, cut),
				  PFHT1050 =normalize(events.HLT.PFHT1050,cut),
				  Photon200=normalize(events.HLT.Photon200,cut),
#				  weight=weight
			)
			#if systematic is None:
			#	variables = {
			#		'phopt':				  leading_pho.pt,
			#		'phophi':				 leading_pho.phi,
			#		'phoeta':				 leading_pho.eta,
			#	}
			#	for variable in output:
			#		if variable not in variables:
			#			continue
			#		normalized_variable = {variable: normalize(variables[variable],cut)}
			#		output[variable].fill(
			#			region=region,
			#			**normalized_variable,
			##			weight=weight,
			#		)


		if shift_name is None:
			systematics = [None] + list(weights.variations)
		else:
			systematics = [shift_name]
			
		for region in regions:
			if region not in selected_regions: continue

			for systematic in systematics:
				if isData and systematic is not None:
					continue
				fill(region, systematic)


		scale = 1
		if self._xsec[dataset]!= -1: 
			scale = 1.0#self._lumi*self._xsec[dataset]

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
