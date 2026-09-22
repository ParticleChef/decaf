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
		'2023pre' : 17.794,
		'2023post': 9.451,
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
			'sr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','JetMET',
					'sig'),
			'wmcr' :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','JetMET',
					'sig'),
			'wecr' :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','EGamma',
					'sig'),
			'tmcr' :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','JetMET',
					'sig'),
			'tecr' :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','EGamma',
					'sig'),
			'zmcr' :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','JetMET',
					'sig'),
			'zecr' :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','EGamma',
					'sig'),
			'gcr'  :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','ST','WW','ZZ','WZ','QCD','EGamma',
					'sig'),
		}
		
		self._TvsQCDwp = {
			'2022pre': 0.33,
			'2022post': 0.33,
			'2023pre': 0.33,
			'2023post': 0.33,
		}
		self._met_triggers = { ## Name is same 22/23 
			'2022pre': [
				'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
				'PFMETNoMu120_PFMHTNoMu120_IDTight'
			],
		}
		self._electron_triggers = { ## Name is same 22/23
			'2022pre':
				[
				'Ele30_WPTight_Gsf',
				'Photon200',
			],
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
			'TvsQCD': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(15,0,1, name='TvsQCD', label='TvsQCD'),
				storage=hist.storage.Weight(),
			),
			'TvsQCD_cat': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(15,0,1, name='TvsQCD', label='TvsQCD'),
				hist.axis.Integer(0,10, name='n_qinfj', label='Number of q from top in AK15 Jet'),
				hist.axis.Integer(0,10, name='n_binfj', label='Number of b flavor in AK15 Jet'),
				storage=hist.storage.Weight(),
			),
			'probT': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(15,0,1, name='probT', label='probT'),
				storage=hist.storage.Weight(),
			),
			'probQCD': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(15,0,1, name='probQCD', label='probQCD'),
				storage=hist.storage.Weight(),
			),
			'fjpt_cat': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Variable([250, 325, 400, 600, 2000], name='fjpt', label='AK15 Leading Jet $p_{T}$', flow=True),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				hist.axis.Integer(0,10, name='n_qinfj', label='Number of q from top in AK15 Jet'),
				hist.axis.Integer(0,10, name='n_binfj', label='Number of b flavor in AK15 Jet'),
				storage=hist.storage.Weight(),
			),
			'fjpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				#hist.axis.Regular(38,250,1200, name='fjpt', label='AK15 Leading Jet $p_{T}$'),
				hist.axis.Variable([250, 325, 400, 600, 2000], name='fjpt', label='AK15 Leading Jet $p_{T}$', flow=True),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
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

		get_met_xy_correction	= self._corrections['get_met_xy_correction']
		get_jec_correction = self._corrections['get_jec_correction']
		get_fjec_correction = self._corrections['get_fjec_correction']
		
		isLooseMuon	 = self._ids['isLooseMuon']	 
		isTightMuon	 = self._ids['isTightMuon']	 
		isLooseElectron = self._ids['isLooseElectron'] 
		isTightElectron = self._ids['isTightElectron'] 
		isLoosePhoton   = self._ids['isLoosePhoton']   
		isTightPhoton   = self._ids['isTightPhoton']   
		isGoodAK4	   = self._ids['isGoodAK4']
		isGoodAK15	= self._ids['isGoodAK15']	
		isJetVeto	= self._ids['isJetVeto']	

		PNetUParTWPs = self._common['btagWPs']['PNetUParT'][self._year]

		###
		#Initialize global quantities (MET ecc.)
		###

		npv = events.PV.npvsGood
		run = events.run
		met = events.MET
		met['pt'] , met['phi'] = get_met_xy_correction(self._year, 'MET', isData,  met.pt, met.phi, npv)

		###
		#Initialize physics objects
		###

		mu = events.Muon
		mu['isloose'] = isLooseMuon(mu,self._year)
		mu['istight'] = isTightMuon(mu,self._year)
		mu_loose=mu[mu.isloose]
		mu_tight=mu[mu.istight]
		mu_ntot = ak.num(mu, axis=1)
		mu_nloose = ak.num(mu_loose, axis=1)
		mu_ntight = ak.num(mu_tight, axis=1)
		
		mu['T'] = ak.zip(
			{
				"r": mu.pt,
				"phi": mu.phi,
			},
			with_name="PolarTwoVector",
			behavior=vector.behavior,
		)
		# define leading mu
		leading_mu = ak.firsts(mu_tight)
		
		# define second mu for Z->mumu
		second_mu = ak.pad_none(mu_tight, target=2)[:,1]
		
		#dimu = ak.cartesian({"mu1": leading_mu, "mu2": second_mu}, nested=True)
		dimu = leading_mu + second_mu
		
		#dimu_mass = (dimu.mu1 + dimu.mu2).mass
		dimu_mass = dimu.mass

		e = events.Electron
		e['isloose'] = isLooseElectron(e,self._year)
		e['istight'] = isTightElectron(e,self._year)

		e['T'] = ak.zip(
			{
				"r": e.pt,
				"phi": e.phi,
			},
			with_name="PolarTwoVector",
			behavior=vector.behavior,
		)

		e_loose = e[e.isloose]
		e_tight = e[e.istight]
		e_ntot = ak.num(e, axis=1)
		e_nloose = ak.num(e_loose, axis=1)
		e_ntight = ak.num(e_tight, axis=1)
		leading_e = ak.firsts(e_tight)

		# define second ele for Z->eleele
		second_e = ak.pad_none(e_tight, target=2)[:,1]
		
		#diele = ak.cartesian({"ele1": leading_ele, "ele2": second_ele}, nested=True)
		diele = leading_e + second_e
		
		#diele_mass = (diele.ele1 + diele.ele2).mass
		diele_mass = diele.mass

		pho = events.Photon
		pho['isloose'] = isLoosePhoton(pho,self._year)
		pho['istight'] = isTightPhoton(pho,self._year)

		pho_loose=pho[pho.isloose]
		pho_tight=pho[pho.istight]
		pho_ntot = ak.num(pho, axis=1)
		pho_nloose = ak.num(pho_loose, axis=1)
		pho_ntight = ak.num(pho_tight, axis=1)
		leading_pho = ak.firsts(pho_tight)
		
		fj = events.AK15PuppiJet
		rho_density = events.Rho.fixedGridRhoFastjetAll
		fjec_corr = get_fjec_correction(self._year, fj.pt, fj.eta, fj.phi, rho_density, fj.area, run, isData)
		fj['pt'] = fj.pt * fjec_corr
		fj['mass'] = fj.mass * fjec_corr
		fj['vec'] = ak.zip(
			{
				"pt": fj.pt,
				"eta": fj.eta,
				"phi": fj.phi,
				"mass": fj.mass,
			},
			with_name="PtEtaPhiMLorentzVector",
			behavior=vector.behavior,
		)
		fj['isclean'] = (
			ak.all(fj.metric_table(mu_loose) > 1.5, axis=2)
			& ak.all(fj.metric_table(e_loose) > 1.5, axis=2)
			& ak.all(fj.metric_table(pho_loose) > 1.5, axis=2)
		)
		fj['isgood'] = isGoodAK15(fj)
		fj['T'] = ak.zip(
			{
				"r": fj.pt,
				"phi": fj.phi,
			},
			with_name="PolarTwoVector",
			behavior=vector.behavior,
		)
		probQCD=fj.ParT_probQCDbb+fj.ParT_probQCDcc+fj.ParT_probQCDb+fj.ParT_probQCDc+fj.ParT_probQCDothers
		probT=fj.ParT_probTopbWqq+fj.ParT_probTopbWcs
		fj['TvsQCD'] = probT/(probT+probQCD)
		fj['probQCD'] = probQCD
		fj['probT'] = probT
		fj_good = fj[fj.isgood]
		fj_clean = fj_good[fj_good.isclean]
		fj_ntot = ak.num(fj, axis=1)
		fj_ngood = ak.num(fj_good, axis=1)
		fj_nclean = ak.num(fj_clean, axis=1)
		leading_fj = ak.firsts(fj_clean)

		j = events.Jet
		jec_corr = get_jec_correction(self._year, j.pt, j.eta, j.phi, rho_density, j.area, run, isData)
		j['pt'] = j.pt * jec_corr
		j['mass'] = j.mass * jec_corr
		j['isveto'] = isJetVeto(j, self._year)
		j_veto = j[j.isveto]
		n_j_veto = ak.num(j_veto, axis=1)
		j['T'] = ak.zip(
			{
				"r": j.pt,
				"phi": j.phi,
			},
			with_name="PolarTwoVector",
			behavior=vector.behavior,
		)

		j['isgood'] = isGoodAK4(j, self._year)
		j['isclean'] = (
			ak.all(j.metric_table(mu_loose) > 0.4, axis=2)
			& ak.all(j.metric_table(e_loose) > 0.4, axis=2)
		)
		j['isiso'] = ak.all(j.metric_table(leading_fj) > 1.5, axis=2)
		j['isbtagvL'] = (j.btagPNetB>PNetUParTWPs['loose'])

		j_good = j[j.isgood]
		j_clean = j_good[j_good.isclean]
		j_iso = j_clean[j_clean.isiso]
		j_btagvL = j_iso[j_iso.isbtagvL]
		leading_j = ak.firsts(j_clean)

		j_ntot=ak.num(j, axis=1)
		j_ngood=ak.num(j_good, axis=1)
		j_nclean=ak.num(j_clean, axis=1)
		j_nbtagvL=ak.num(j_btagvL, axis=1)



		if not isData:
			j['isbflav'] = j.hadronFlavour == 5
			j_bflav = j[j.isbflav]
			#fj['bisinfj'] = ak.all(fj.vec.metric_table(j_bflav)<1.5, axis=2)
			binfj = leading_fj.vec.delta_r(j_bflav)
			n_binfj = ak.sum(binfj < 1.5, axis=1)
			#fj_bflav = fj[fj.bisinfj]
			#nfj_bflav = ak.num(fj_bflav, axis=1)
			
			gen = events.GenPart

			yourmotheris = gen.distinctParent.distinctParent.pdgId
			gen['isq'] = (abs(gen.pdgId) >= 1) & (abs(gen.pdgId) <= 5) & gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			p1 = ak.fill_none(gen.distinctParent.pdgId, 0)
			p2 = ak.fill_none(ak.fill_none(gen.distinctParent.distinctParent.pdgId, 0), 0)
			gen['isqfromt'] = gen.isq & (abs(p1) == 24) & (abs(p2) == 6)
			gen['isbfromt'] = gen.isq & (abs(p1) == 6)
			q_from_t = gen[gen.isqfromt]
			b_from_t = gen[gen.isbfromt]
			qinfj = leading_fj.vec.delta_r(q_from_t)
			n_qinfj = ak.sum(qinfj < 1.5, axis=1)
			#fj['qisinfj'] = ak.all(fj.vec.metric_table(q_from_t)<1.5, axis=2)
			#q_from_t_fj = fj[fj.qisinfj]
			#nq_from_t_fj = ak.num(q_from_t_fj, axis=1)

			#gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			#gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			#genTops = gen[gen.isTop]
			#genW = gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))]
			#gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			GenIsoPho = events.GenIsolatedPhoton
			nGenIsoPho = ak.num(GenIsoPho, axis=1)
		
		###
		# Selections
		###

		lumimask = np.ones(len(events), dtype='bool')
		if isData:
			lumimask = AnalysisProcessor.lumiMasks[self._year](events.run, events.luminosityBlock)
		selection.add('lumimask', lumimask)

		met_filters =  np.ones(len(events), dtype='bool')
		#for flag in AnalysisProcessor.met_filters[self._year]:
		for flag in AnalysisProcessor.met_filters['2022pre']:
			met_filters = met_filters & events.Flag[flag]
		selection.add('met_filters',met_filters)

		triggers = np.zeros(len(events), dtype='bool')
		for path in self._met_triggers['2022pre']:
			if not hasattr(events.HLT, path): continue
			triggers = triggers | events.HLT[path]
		selection.add('met_triggers', triggers)

		triggers = np.zeros(len(events), dtype='bool')
		for path in self._electron_triggers['2022pre']:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('single_electron_triggers', ak.to_numpy(triggers))

		triggers = np.zeros(len(events), dtype='bool')
		for path in self._photon_triggers['2022pre']:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('single_photon_triggers', ak.to_numpy(triggers))

		if ('WtoLNu-2Jets' in dataset) & ('PT' in dataset):
			remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt >120) ## W
			selection.add("exclude_wjets_greater_120", ak.to_numpy(ak.all(remove_overlap, axis=1)))
		else:
			selection.add("exclude_wjets_greater_120", np.full(len(events), True))

		if ('WtoLNu-2Jets' in dataset) & (not ('PT' in dataset)):
			remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt <= 120) ## w
			selection.add("exclude_wjets_less_120", ak.to_numpy(ak.all(remove_overlap, axis=1)))
		else:
			selection.add("exclude_wjets_less_120", np.full(len(events), True))

		if 'QCD_PT' in dataset:
			selection.add("QCD_NoGenIsoPho", (nGenIsoPho==0))
		else:
			selection.add("QCD_NoGenIsoPho", np.full(len(events), True))
		
		if 'GJ_PTG' in dataset:
			selection.add("GJet_GenIsoPho", (nGenIsoPho>=1))
		else:
			selection.add("GJet_GenIsoPho", np.full(len(events), True))
		

		
		selection.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(pho_nloose==0))
		selection.add('isoneM',  (e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(pho_nloose==0))
		selection.add('isoneE',  (e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(pho_nloose==0))
		selection.add('isoneG',  (e_nloose==0)&(mu_nloose==0)&(pho_nloose==1)&(pho_ntight==1))
		selection.add('istwoM',  (e_nloose==0)&(mu_nloose==2)&(pho_nloose==0))
		selection.add('istwoE',  (e_nloose==2)&(mu_nloose==0)&(pho_nloose==0))
		selection.add('one_ak4',  (j_nclean>0))
		selection.add('one_ak15', (fj_nclean>0))
		selection.add('leading_fj250', (leading_fj.pt>250))

		selection.add('noextrab', (j_nbtagvL==0))
		selection.add('extrab',   (j_nbtagvL>0))
		selection.add('jetveto', (n_j_veto==0))

		selection.add('met120', (met.pt<120))
		selection.add('met150', (met.pt>150))
		selection.add('diele_mass', (diele_mass>60)&(diele_mass<120))
		selection.add('dimu_mass',  (dimu_mass>60)&(dimu_mass<120))



		regions = {
			'sr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'iszeroL',
					'noextrab',
			],
			'wmcr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'isoneM',
					'met150',
					'noextrab',
			],
			'wecr': [
					'lumimask',
					'met_filters', 'single_electron_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'isoneE',
					'met150',
					'noextrab',
			],
			'tmcr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'isoneM',
					'met150',
					'extrab',
			],
			'tecr': [
					'lumimask',
					'met_filters', 'single_electron_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'isoneE',
					'met150',
					'extrab',
			],
			'zmcr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'istwoM',
					'met120',
					'dimu_mass',
			],
			'zecr': [
					'lumimask',
					'met_filters', 'single_electron_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'istwoE',
					'met120',
					'diele_mass',
			],
			'gcr': [
					'lumimask',
					'met_filters', 'single_photon_triggers',
					'exclude_wjets_greater_120', 'exclude_wjets_less_120',
					'QCD_NoGenIsoPho', 'GJet_GenIsoPho',
					'jetveto',
					'one_ak15',
					'leading_fj250',
					'isoneG',
					'noextrab',
			],
		}

		def normalize(val, cut, fillvalue=np.nan):
			# fillvalue must stay integer (e.g. -1) for values feeding an Integer axis,
			# otherwise ak.fill_none(..., np.nan) upcasts the whole array to float64
			# and boost-histogram rejects it with "Only integer arrays supported
			# when targeting integer axes".
			if cut is None:
				ar = ak.to_numpy(ak.fill_none(val, fillvalue))
				return ar
			else:
				ar = ak.to_numpy(ak.fill_none(val[cut], fillvalue))
				return ar
				
		def fill(region, systematic):
			cut = selection.all(*regions[region])
			if systematic in weights.variations:
				weight = weights.weight(modifier=systematic)[cut]
			else:
				weight = weights.weight()[cut]
			if  not isData:
				output['fjpt_cat'].fill(
					  region=region,
					  fjpt=normalize(leading_fj.pt, cut),
					  TvsQCD=normalize(leading_fj.TvsQCD, cut),
					  n_qinfj=normalize(n_qinfj, cut, fillvalue=-1),
					  n_binfj=normalize(n_binfj, cut, fillvalue=-1),
					  weight=weight
				)
				output['TvsQCD_cat'].fill(
					  region=region,
					  TvsQCD=normalize(leading_fj.TvsQCD, cut),
					  n_qinfj=normalize(n_qinfj, cut, fillvalue=-1),
					  n_binfj=normalize(n_binfj, cut, fillvalue=-1),
					  weight=weight
				)
			output['fjpt'].fill(
				  region=region,
				  fjpt=normalize(leading_fj.pt, cut),
				  TvsQCD=normalize(leading_fj.TvsQCD, cut),
				  weight=weight
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
			scale = self._lumi*self._xsec[dataset]
			#scale = 1.0#self._lumi*self._xsec[dataset]

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
