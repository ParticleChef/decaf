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
			'basic' :('WtoLNu-2Jets','DYto2L-2Jets','Muon_Run2022C','Muon_Run2022D'),
			'muon' :('WtoLNu-2Jets','Muon_Run2022C','Muon_Run2022D'),
			'dimu' :('DYto2L-2Jets','Muon_Run2022C','Muon_Run2022D'),
		}
		

		self._met_triggers = { ## NOT UPDATED for Run3
			'2022pre': [
				'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
				'PFMETNoMu120_PFMHTNoMu120_IDTight'
			]
		}
		self._met_ref_triggers = {
			'2022pre':
				[
				'IsoMu24'
			]
		}
		self._corrections = corrections
		self._ids = ids
		self._common = common

		self.make_output = lambda: {
			'sumw': 0.,
			'template': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='mupt', label='Muon Pt'),
				hist.axis.Variable([-0.5,0.5,1.5], name='IsoMu24', label='Reference trigger'),#growth=False),
				hist.axis.Variable([-0.5,0.5,1.5], name='PFMETNoMu120',label='Photon trigger'),#growth=False),
				#hist.axis.Variable([40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300], name='fjmass', label=r'AK15 Jet $m_{sd}$'),
				storage=hist.storage.Weight(),
			),
			'mupt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(24,0,1200, name='mupt', label='Leading Muon Pt'),
				hist.axis.Integer(0,2, name='IsoMu', label='Reference trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120', label='PFMETNoMu120 trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120PFHT60', label='PFMETNoMu120PFHT60 trigger'),
				storage=hist.storage.Weight(),
			),
			'mueta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='mueta', label='Leading Muon Eta'),
				hist.axis.Integer(0,2, name='IsoMu', label='Reference trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120', label='PFMETNoMu120 trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120PFHT60', label='PFMETNoMu120PFHT60 trigger'),
				storage=hist.storage.Weight(),
			),
			'muphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='muphi', label='Leading Muon Phi'),
				hist.axis.Integer(0,2, name='IsoMu', label='Reference trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120', label='PFMETNoMu120 trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120PFHT60', label='PFMETNoMu120PFHT60 trigger'),
				storage=hist.storage.Weight(),
			),
			'recoilpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(24,0,1200, name='recoilpt', label='Leading Recoil Pt'),
				hist.axis.Integer(0,2, name='IsoMu', label='Reference trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120', label='PFMETNoMu120 trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120PFHT60', label='PFMETNoMu120PFHT60 trigger'),
				storage=hist.storage.Weight(),
			),
			'recoilphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='recoilphi', label='Leading Recoil Phi'),
				hist.axis.Integer(0,2, name='IsoMu', label='Reference trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120', label='PFMETNoMu120 trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120PFHT60', label='PFMETNoMu120PFHT60 trigger'),
				storage=hist.storage.Weight(),
			),
			'metpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(24,0,1200, name='metpt', label='MET Pt'),
				hist.axis.Integer(0,2, name='IsoMu', label='Reference trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120', label='PFMETNoMu120 trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120PFHT60', label='PFMETNoMu120PFHT60 trigger'),
				storage=hist.storage.Weight(),
			),
			'metphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='metphi', label='MET Phi'),
				hist.axis.Integer(0,2, name='IsoMu', label='Reference trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120', label='PFMETNoMu120 trigger'),
				hist.axis.Integer(0,2, name='PFMETNoMu120PFHT60', label='PFMETNoMu120PFHT60 trigger'),
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

		get_met_xy_correction	= self._corrections['get_met_xy_correction']

		###
		#Initialize global quantities (MET ecc.)
		###

		npv = events.PV.npvsGood
		met = events.MET
		met['pt'] , met['phi'] = get_met_xy_correction(self._year, 'MET', isData,  met.pt, met.phi, npv)

		###
		#Initialize physics objects
		###

		mu = events.Muon
		mu['isloose'] = isLooseMuon(mu,self._year)
		mu['istight'] = isTightMuon(mu,self._year)
		mu['T'] = ak.zip(
			{
				"r": mu.pt,
				"phi": mu.phi,
			},
			with_name="PolarTwoVector",
			behavior=vector.behavior,
		)

		mu['ispt30'] = (ak.firsts(mu).pt > 30)
		mu_pt30 = mu[mu.ispt30]
		mu_loose = mu[mu.isloose]
		mu_tight = mu[mu.istight]
		mu_npt30 = ak.num(mu_pt30, axis=1)
		mu_nloose = ak.num(mu_loose, axis=1)
		mu_ntight = ak.num(mu_tight, axis=1)
		
		leading_mu = ak.firsts(mu)
		second_mu = ak.pad_none(mu, target=2)[:,1]
		
		dimu = leading_mu + second_mu
		dimu_mass = dimu.mass


		u = {
			'basic'  : met+leading_mu.T,
			'muon'  : met+leading_mu.T,
			'dimu'  : met+leading_mu.T+second_mu.T,
		}

		
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
		selection.add('met_triggers', triggers)

		## Reference trigger
		triggers = np.zeros(len(events), dtype='bool')
		for path in self._met_ref_triggers[self._year]:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('reference_triggers', ak.to_numpy(triggers))


		

		selection.add('met100', (met.pt > 100))
		selection.add('met120', (met.pt < 120))
		selection.add('muon_selection', ((mu_ntight==1)&(mu_nloose==1)))
		selection.add('dimuon_selection', (mu_npt30==2))
		selection.add('mu1pt30', (ak.firsts(mu).pt > 30))
		selection.add('dimu_mass',(dimu_mass>60)&(dimu_mass<120))
		selection.add('recoil_muon', (u['muon'].pt > 70))
		selection.add('recoil_dimu', (u['dimu'].pt > 70))


		regions = {
			'basic': [
					'lumimask',
					'met_filters',
			],
			'muon': [
					'lumimask',
					'met_filters',
					'recoil_muon',
					'muon_selection',
					'met100',
			],
			'dimu': [
					'lumimask',
					'met_filters',
					'recoil_dimu',
					'dimuon_selection',
					'met120',
					'mu1pt30',
					'dimu_mass',
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
				  mupt=normalize(leading_mu.pt, cut),
				  IsoMu24 =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
			)
			output['mupt'].fill(
				  region=region,
				  mupt=normalize(leading_mu.pt, cut),
				  IsoMu =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
				  PFMETNoMu120PFHT60=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,cut),
			)
			output['muphi'].fill(
				  region=region,
				  muphi=normalize(leading_mu.phi, cut),
				  IsoMu =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
				  PFMETNoMu120PFHT60=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,cut),
			)
			output['mueta'].fill(
				  region=region,
				  mueta=normalize(leading_mu.eta, cut),
				  IsoMu =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
				  PFMETNoMu120PFHT60=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,cut),
			)
			output['recoilpt'].fill(
				  region=region,
				  recoilpt=normalize(u[region].pt, cut),
				  IsoMu =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
				  PFMETNoMu120PFHT60=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,cut),
			)
			output['recoilphi'].fill(
				  region=region,
				  recoilphi=normalize(u[region].phi, cut),
				  IsoMu =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
				  PFMETNoMu120PFHT60=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,cut),
			)
			output['metpt'].fill(
				  region=region,
				  metpt=normalize(met.pt, cut),
				  IsoMu =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
				  PFMETNoMu120PFHT60=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,cut),
			)
			output['metphi'].fill(
				  region=region,
				  metphi=normalize(met.phi, cut),
				  IsoMu =normalize(events.HLT.IsoMu24,cut),
				  PFMETNoMu120=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight,cut),
				  PFMETNoMu120PFHT60=normalize(events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,cut),
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
