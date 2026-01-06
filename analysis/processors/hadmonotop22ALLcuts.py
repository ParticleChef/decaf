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
		'2022pre' : 7.980,
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
			'sr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','JetMET'),
			'wmcr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','JetMET'),
			'wecr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','EGamma'),
			'tmcr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','JetMET'),
			'tecr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','EGamma'),
			'zmcr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','JetMET'),
			'zecr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','EGamma'),
			'gcr'   :('TT','WtoLNu-2Jets','DYto2L-2Jets','Zto2Nu-2Jets','GJ','TBbar','TQ','Tbar','TW','WW','ZZ','WZ','QCD','EGamma'),
		}
		

		self._met_triggers = { 
			'2022pre': [
				'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
				'PFMETNoMu120_PFMHTNoMu120_IDTight'
			]
		}
		self._electron_triggers = {
			'2022pre':
				[
				'Ele30_WPTight_Gsf',
				'Photon200',
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
			'cutflow': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.StrCategory([], name='cutname', growth=True),
				hist.axis.Variable([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], name='cutflow', label='cut'),
				storage=hist.storage.Weight(),
			),
			'template': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.StrCategory([], name='systematic', growth=True),
				#hist.axis.Variable([40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300], name='fjmass', label=r'AK15 Jet $m_{sd}$'),
				storage=hist.storage.Weight(),
			),
			'TvsQCD': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(15,0,1, name='TvsQCD', label='TvsQCD'),
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
			'mindphirecoil': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(30,0,3.5, name='mindphirecoil', label='Min |dPhi(Recoil,AK4s)|'),
				storage=hist.storage.Weight(),
			),
			'minDphirecoil': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(30,0,3.5, name='minDphirecoil', label='Min |dPhi(Recoil, leading AK15s)|'),
				storage=hist.storage.Weight(),
			),
			'ut': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(65,350,1000, name='ut', label='U_{T}'),
				storage=hist.storage.Weight(),
			),
			'uphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='uphi', label='U_{phi}'),
				storage=hist.storage.Weight(),
			),
			'mupt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='mupt', label='Leading Muon Pt'),
				storage=hist.storage.Weight(),
			),
			'mueta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='mueta', label='Leading Muon Eta'),
				storage=hist.storage.Weight(),
			),
			'muphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='muphi', label='Leading Muon Phi'),
				storage=hist.storage.Weight(),
			),
			'elept': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='elept', label='Leading Electron Pt'),
				storage=hist.storage.Weight(),
			),
			'eleeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='eleeta', label='Leading Electron Eta'),
				storage=hist.storage.Weight(),
			),
			'elephi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='elephi', label='Leading Electron Phi'),
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
			'dielept': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='dielept', label='DiElectron Pt'),
				storage=hist.storage.Weight(),
			),
			'dieleeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='dieleeta', label='DiElectron Eta'),
				storage=hist.storage.Weight(),
			),
			'dielephi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='dielephi', label='DiElectron Phi'),
				storage=hist.storage.Weight(),
			),
			'dimupt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='dimupt', label='DiMuon Pt'),
				storage=hist.storage.Weight(),
			),
			'dimueta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='dimueta', label='DiMuon Eta'),
				storage=hist.storage.Weight(),
			),
			'dimuphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='dimuphi', label='DiMuon Phi'),
				storage=hist.storage.Weight(),
			),
			'jbtagLpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='jbtagLpt', label='AK4 Leading BTag Jet $p_{T}$'),
				storage=hist.storage.Weight(),
			),
			'jbtagLeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='jbtagLeta', label='AK4 Leading BTag Jet Eta'),
				storage=hist.storage.Weight(),
			),
			'jbtagLphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='jbtagLphi', label='AK4 Leading BTag Jet Phi'),
				storage=hist.storage.Weight(),
			),
			'jpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='jpt', label='AK4 Leading Jet $p_{T}$'),
				storage=hist.storage.Weight(),
			),
			'jeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='jeta', label='AK4 Leading Jet Eta'),
				storage=hist.storage.Weight(),
			),
			'jphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='jphi', label='AK4 Leading Jet Phi'),
				storage=hist.storage.Weight(),
			),
			'fjpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='fjpt', label='AK15 Leading Jet Pt'),
				storage=hist.storage.Weight(),
			),
			'fjeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='fjeta', label='AK15 Leading Jet Eta'),
				storage=hist.storage.Weight(),
			),
			'fjphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='fjphi', label='AK15 Leading Jet Phi'),
				storage=hist.storage.Weight(),
			),
			'mT': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,1200, name='mT', label='mT'),
				storage=hist.storage.Weight(),
			),
			'met': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(30,0,1200, name='met', label='MET $p_{T}$ [GeV]'),
				storage=hist.storage.Weight(),
			),
			'metphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='metphi', label='MET $\phi$'),
				#hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'njets': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='njets', label='AK4 Number of Jets'),
				storage=hist.storage.Weight(),
			),
			'njbtagL': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='njbtagL', label='AK4 Number of BTag Jets'),
				storage=hist.storage.Weight(),
			),
			'nfjets': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='nfjets', label='AK15 Number of Jets'),
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
	

		#jet_factory			  = self._corrections['jet_factory']
		#fatjet_factory		   = self._corrections['fatjet_factory']
		#subjet_factory		   = self._corrections['subjet_factory']
		#met_factory			  = self._corrections['met_factory']

		#
		#jec_cache = cachetools.Cache(np.inf)
	
		#nojer = "NOJER" if self._skipJER else ""
		#thekey = f"{self._year}mc{nojer}"

		#def add_jec_variables(jets, event_rho):
		#	jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
		#	jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
		#	jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
		#	jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
		#	return jets
		#
		#def add_jec_variables_fj(jets, event_rho):
		#	jets["pt_raw"] = (1 - jets.rawFactor)*jets.Jet_pt
		#	jets["mass_raw"] = (1 - jets.rawFactor)*jets.Jet_mass
		#	jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
		#	jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.Jet_pt)[0]
		#	return jets

		#jets = jet_factory[thekey].build(add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll), jec_cache)
		#fatjets = fatjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiJet, events.fixedGridRhoFastjetAll), jec_cache)
		#subjets = subjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiSubJet, events.fixedGridRhoFastjetAll), jec_cache)
		#met = met_factory.build(events.MET, jets, {})

		#shifts = [({"Jet": jets, "AK15PFPuppiSubJet": subjets, "AK15PFPuppiJet": fatjets, "MET": met}, None)]
		#if self._systematics:
		#	shifts.extend([
		#		({"Jet": jets.JES_jes.up, "AK15PFPuppiSubJet": subjets.JES_jes.up, "AK15PFPuppiSubJet": fatjets.JES_jes.up, "MET": met.JES_jes.up}, "JESUp"),
		#		({"Jet": jets.JES_jes.down, "AK15PFPuppiSubJet": subjets.JES_jes.down, "AK15PFPuppiSubJet": fatjets.JES_jes.down, "MET": met.JES_jes.down}, "JESDown"),
		#		({"Jet": jets, "AK15PFPuppiSubJet": subjets, "AK15PFPuppiSubJet": fatjets, "MET": met.MET_UnclusteredEnergy.up}, "UESUp"),
		#		({"Jet": jets, "AK15PFPuppiSubJet": subjets, "AK15PFPuppiSubJet": fatjets, "MET": met.MET_UnclusteredEnergy.down}, "UESDown"),
		#	])
		#	if not self._skipJER:
		#		shifts.extend([
		#			({"Jet": jets.JER.up, "AK15PFPuppiSubJet": subjets.JER.up, "AK15PFPuppiSubJet": fatjets.JER.up, "MET": met.JER.up}, "JERUp"),
		#			({"Jet": jets.JER.down, "AK15PFPuppiSubJet": subjets.JER.down, "AK15PFPuppiSubJet": fatjets.JER.down, "MET": met.JER.down}, "JERDown"),
		#		])
		#return processor.accumulate(self.process_shift(update(events, collections), name) for collections, name in shifts)

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

		get_pu_weight			= self._corrections['get_pu_weight']	
		get_mu_highpt_id_sf	   = self._corrections['get_mu_highpt_id_sf']
		get_mu_hlt_sf		 = self._corrections['get_mu_hlt_sf']
		get_mu_sf				= self._corrections['get_mu_sf']
		get_mu_loose_id_sf	   = self._corrections['get_mu_loose_id_sf']
		get_mu_tight_id_sf	   = self._corrections['get_mu_tight_id_sf']
		get_mu_loose_iso_sf	  = self._corrections['get_mu_loose_iso_sf']
		get_mu_tight_iso_sf	  = self._corrections['get_mu_tight_iso_sf']
		#get_met_trig_weight	  = self._corrections['get_met_trig_weight']
		#get_ele_loose_id_sf	  = self._corrections['get_ele_loose_id_sf']
		#get_ele_tight_id_sf	  = self._corrections['get_ele_tight_id_sf']
		#get_ele_trig_weight	  = self._corrections['get_ele_trig_weight']
		#get_ele_reco_sf_below20  = self._corrections['get_ele_reco_sf_below20']
		#get_ele_reco_sf_above20  = self._corrections['get_ele_reco_sf_above20']
		#get_mu_rochester_sf	  = self._corrections['get_mu_rochester_sf'][self._year]
		#get_met_xy_correction	= self._corrections['get_met_xy_correction']
		#get_nlo_ewk_weight	   = self._corrections['get_nlo_ewk_weight']	
		#get_nnlo_nlo_weight	  = self._corrections['get_nnlo_nlo_weight'][self._year]
		#get_btag_weight	  = self._corrections['get_btag_weight']
		#get_ttbar_weight	 = self._corrections['get_ttbar_weight']
		
		isLooseMuon	 = self._ids['isLooseMuon']	 
		isTightMuon	 = self._ids['isTightMuon']	 
		isLooseElectron = self._ids['isLooseElectron'] 
		isTightElectron = self._ids['isTightElectron'] 
		isLoosePhoton   = self._ids['isLoosePhoton']   
		isTightPhoton   = self._ids['isTightPhoton']   
		isGoodAK4	   = self._ids['isGoodAK4']	   
		isGoodAK15	= self._ids['isGoodAK15']	
		
		PNetUParTWPs = self._common['btagWPs']['PNetUParT'][self._year]


		###
		#Initialize global quantities (MET ecc.)
		###

		npv = events.PV.npvsGood
		run = events.run
		met = events.MET
		#met['pt'] , met['phi'] = get_met_xy_correction(self._year, npv, run, met.pt, met.phi, isData)

		###
		#Initialize physics objects
		###

		mu = events.Muon
		mu['isloose'] = isLooseMuon(mu,self._year)
		mu['id_sf'] = ak.where(
			mu.isloose, 
			get_mu_loose_id_sf(self._year, abs(mu.eta), mu.pt), 
			ak.ones_like(mu.pt)
		)
		#mu['iso_sf'] = ak.where(
		#	mu.isloose, 
		#	get_mu_loose_iso_sf(self._year, abs(mu.eta), mu.pt), 
		#	ak.ones_like(mu.pt)
		#)
		mu['istight'] = isTightMuon(mu,self._year)
		#mu['id_sf'] = ak.where(
		#	mu.istight, 
		#	get_mu_tight_id_sf(self._year, abs(mu.eta), mu.pt), 
		#	mu.id_sf
		#)
		#mu['iso_sf'] = ak.where(
		#	mu.istight, 
		#	get_mu_tight_iso_sf(self._year, abs(mu.eta), mu.pt), 
		#	ak.ones_like(mu.pt) #mu.iso_sf
		#)
		mu['T'] = ak.zip(
			{
				"r": mu.pt,
				"phi": mu.phi,
			},
			with_name="PolarTwoVector",
			behavior=vector.behavior,
		)
		mu_loose=mu[mu.isloose]
		mu_tight=mu[mu.istight]
		mu_ntot = ak.num(mu, axis=1)
		mu_nloose = ak.num(mu_loose, axis=1)
		mu_ntight = ak.num(mu_tight, axis=1)
		
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
		pho['T'] = ak.zip(
			{
				"r": pho.pt,
				"phi": pho.phi,
			},
			with_name="PolarTwoVector",
			behavior=vector.behavior,
		)

		pho_loose=pho[pho.isloose]
		pho_tight=pho[pho.istight]
		pho_ntot = ak.num(pho, axis=1)
		pho_nloose = ak.num(pho_loose, axis=1)
		pho_ntight = ak.num(pho_tight, axis=1)
		leading_pho = ak.firsts(pho_tight)

		fj = events.AK15PuppiJet
		#fj['msd_corr'] = get_msd_corr(fj)
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

		###
		# Calculate recoil and transverse mass
		###

		u = {
			'sr'	: met,			   
			'wmcr'  : met+leading_mu.T,
			'wecr'  : met+leading_e.T,
			'tmcr'  : met+leading_mu.T,
			'tecr'  : met+leading_e.T,
			'zmcr'  : met+leading_mu.T+second_mu.T,
			'zecr'  : met+leading_e.T+second_e.T,
			'gcr'   : met+leading_pho.T 
		}

		Dphi = {
			'sr'	: u[region].delta_phi(leading_fj.T),
			'wmcr'  : u[region].delta_phi(leading_fj.T),
			'wecr'  : u[region].delta_phi(leading_fj.T),
			'tmcr'  : u[region].delta_phi(leading_fj.T),
			'tecr'  : u[region].delta_phi(leading_fj.T),
			'zmcr'  : u[region].delta_phi(leading_fj.T),
			'zecr'  : u[region].delta_phi(leading_fj.T),
			'gcr'   : u[region].delta_phi(leading_fj.T)
		}

		dphi = {
			'sr'	: u['sr'].delta_phi(j_clean.T),
			'wmcr'  : u['sr'].delta_phi(j_clean.T),
			'wecr'  : u['sr'].delta_phi(j_clean.T),
			'tmcr'  : u['sr'].delta_phi(j_clean.T),
			'tecr'  : u['sr'].delta_phi(j_clean.T),
			'zmcr'  : u['sr'].delta_phi(j_clean.T),
			'zecr'  : u['sr'].delta_phi(j_clean.T),
			'gcr'   : u['sr'].delta_phi(j_clean.T)
		}

		mT = {
			#'sr'   : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
			'wmcr' : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
			'wecr' : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
			'tmcr' : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
			'tecr' : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
		}

		###
		#Calculating weights
		###
		if not isData:
			
			gen = events.GenPart

			gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			genTops = gen[gen.isTop]
			#nlo = np.ones(len(events), dtype='float')
			##print('genTops pt ', genTops.pt)
			#if('TTTo' in dataset): 
			#	nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt) * get_ttbar_weight(genTops[:,1].pt))
			#	
			#gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			#gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			#gen['isA'] = (abs(gen.pdgId)==22) & gen.hasFlags(['isPrompt', 'fromHardProcess', 'isLastCopy']) & (gen.status == 1)
			#
			#genWs = gen[gen.isW] 
			#genZs = gen[gen.isZ]
			#genDYs = gen[gen.isZ&(gen.mass>30)]
			#genAs = gen[gen.isA & (gen.pt > 100)]
			#
			#GenIsoPho = events.GenIsolatedPhoton
			#nGenIsoPho = ak.num(GenIsoPho, axis=1)
			
			#nlo_ewk = np.ones(len(events), dtype='float')
			#if('WJets' in dataset): 
			#	#nlo_ewk = get_nlo_ewk_weight['w'](genWs.pt.max())
			#	nlo_ewk = get_nlo_ewk_weight['w'](ak.max(genWs.pt, axis=1))
			#elif('DY' in dataset): 
			#	nlo_ewk = get_nlo_ewk_weight['dy'](ak.max(genDYs.pt, axis=1))
			#elif('Z1Jets' in dataset or 'Z2Jets' in dataset): 
			#	nlo_ewk = get_nlo_ewk_weight['z'](ak.max(genZs.pt, axis=1))
			#elif('G1Jet' in dataset):
			#	nlo_ewk = get_nlo_ewk_weight['a'](ak.max(genAs.pt, axis=1))

			###
			# Calculate PU weight and systematic variations
			###
			#pu = get_pu_weight(self._year, events.Pileup.nTrueInt)
			pu ={
				'sr'  : np.ones(len(events), dtype='float'),
				'wmcr': get_pu_weight(self._year, events.Pileup.nTrueInt),
				'wecr': get_pu_weight(self._year, events.Pileup.nTrueInt),
				'gcr' : get_pu_weight(self._year, events.Pileup.nTrueInt),
			}

			###
			# Trigger efficiency weight
			###
			trig ={
				'sr'  : np.ones(len(events), dtype='float'),
				'wmcr': get_pu_weight(self._year, events.Pileup.nTrueInt),
				'wecr': get_pu_weight(self._year, events.Pileup.nTrueInt),
				'gcr' : get_pu_weight(self._year, events.Pileup.nTrueInt),
			}

			#trig = {
			#	'nocr':   get_met_trig_weight(self._year, met.pt),
			#	'pucr': get_met_trig_weight(self._year, u['wmcr'].r),
			#}

			### 
			# Calculating electron and muon ID weights
			###
			#ids ={
			#	'nocr':  np.ones(len(events), dtype='float'),
			#	'pucr':  np.ones(len(events), dtype='float'),
			#	'allcr': leading_mu.highpt_id#get_mu_highpt_id_sf(self._year, abs(leading_mu.eta), leading_mu.pt),
			#}
		   
			###
			# Reconstruction weights for electrons
			###									   
			#reco = {
			#	'nocr': np.ones(len(events), dtype='float'),
			#	'pucr': np.ones(len(events), dtype='float'),
			#	'allcr': np.ones(len(events), dtype='float'),
			#}

			###
			# Isolation weights for muons
			###
			#isolation = {
			#	'nocr':  np.ones(len(events), dtype='float'),
			#	'pucr':  np.ones(len(events), dtype='float'),
			#	'allcr': leading_mu.iso_sf#get_mu_tight_iso_sf(self._year, abs(leading_mu.eta), leading_mu.pt),
			#}

			###
			# HLT weights for muons
			###
			#hlt= {
			#	'nocr':  np.ones(len(events), dtype='float'),
			#	'pucr':  np.ones(len(events), dtype='float'),
			#	'allcr': leading_mu.hlt_id#get_mu_hlt_sf(self._year, abs(leading_mu.eta), leading_mu.pt),
			#}

			###
			# AK4 b-tagging weights
			###
			#btagSF, \
			#btagSFbc_correlatedUp, \
			#btagSFbc_correlatedDown, \
			#btagSFbc_uncorrelatedUp, \
			#btagSFbc_uncorrelatedDown, \
			#btagSFlight_correlatedUp, \
			#btagSFlight_correlatedDown, \
			#btagSFlight_uncorrelatedUp, \
			#btagSFlight_uncorrelatedDown  = get_btag_weight('deepflav',self._year,'loose').btag_weight(
			#	j_iso.pt,
			#	j_iso.eta,
			#	j_iso.hadronFlavour,
			#	j_iso.isdflvL
			#)

			weights.add('genw',events.genWeight)
			#weights.add('pileup',pu[region])
			#weights.add('ids', ids[region])
			#weights.add('isolation', isolation[region])
			#weights.add('hlt', hlt[region])
			#weights.add('trig', trig[region])
			#weights.add('nlo_ewk',nlo_ewk)
			#weights.add('reco', reco[region])
			#weights.add('btagSF',btagSF)
			#weights.add('btagSFbc_correlated',np.ones(len(events), dtype='float'), btagSFbc_correlatedUp/btagSF, btagSFbc_correlatedDown/btagSF)
			#weights.add('btagSFbc_uncorrelated',np.ones(len(events), dtype='float'), btagSFbc_uncorrelatedUp/btagSF, btagSFbc_uncorrelatedDown/btagSF)
			#weights.add('btagSFlight_correlated',np.ones(len(events), dtype='float'), btagSFlight_correlatedUp/btagSF, btagSFlight_correlatedDown/btagSF)
			#weights.add('btagSFlight_uncorrelated',np.ones(len(events), dtype='float'), btagSFlight_uncorrelatedUp/btagSF, btagSFlight_uncorrelatedDown/btagSF)
			
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

		triggers = np.zeros(len(events), dtype='bool')
		for path in self._electron_triggers[self._year]:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('single_electron_triggers', ak.to_numpy(triggers))

		triggers = np.zeros(len(events), dtype='bool')
		for path in self._photon_triggers[self._year]:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('single_photon_triggers', ak.to_numpy(triggers))

		if ('WtoLNu-2Jets' in dataset) & ('PT' in dataset):
			remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt >120) ## W
			selection.add("exclude_wjets_greater_400", ak.to_numpy(ak.all(remove_overlap, axis=1)))
		else:
			selection.add("exclude_wjets_greater_400", np.full(len(events), True))

		if ('WtoLNu-2Jets' in dataset) & (not ('PT' in dataset)):
			remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt <= 120) ## w
			selection.add("exclude_wjets_less_400", ak.to_numpy(ak.all(remove_overlap, axis=1)))
		else:
			selection.add("exclude_wjets_less_400", np.full(len(events), True))

		

		selection.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(pho_nloose==0))
		selection.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(pho_nloose==0))
		selection.add('isoneE', (e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(pho_nloose==0))
		selection.add('isoneG', (e_nloose==0)&(mu_nloose==0)&(pho_nloose==1)&(pho_ntight==1))
		selection.add('istwoM', (e_nloose==0)&(mu_nloose==2)&(pho_nloose==0))
		selection.add('istwoE', (e_nloose==2)&(mu_nloose==0)&(pho_nloose==0))
		selection.add('one_ak4', (j_nclean>0))
		selection.add('one_ak15', (fj_nclean>0))
		selection.add('leading_fj250', (leading_fj.pt>250))

		selection.add('noextrab', (j_nbtagvL==0))
		selection.add('extrab', (j_nbtagvL>0))

		selection.add('met120',(met.pt<120))
		selection.add('met150',(met.pt>150))
		selection.add('diele_mass',(diele_mass>60)&(diele_mass<120))
		selection.add('dimu_mass',(dimu_mass>60)&(dimu_mass<120))


		regions = {
			'sr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'iszeroL',
					'noextrab',
			],
			'wmcr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'isoneM',
					'met150',
					'noextrab',
			],
			'wecr': [
					'lumimask',
					'met_filters', 'single_electron_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'isoneE',
					'met150',
					'noextrab',
			],
			'tmcr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'isoneM',
					'met150',
					'extrab',
			],
			'tecr': [
					'lumimask',
					'met_filters', 'single_electron_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'isoneE',
					'met150',
					'extrab',
			],
			'zmcr': [
					'lumimask',
					'met_filters', 'met_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'istwoM',
					'met120',
					'dimu_mass',
			],
			'zecr': [
					'lumimask',
					'met_filters', 'single_electron_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'istwoE',
					'met120',
					'diele_mass',
			],
			'gcr': [
					'lumimask',
					'met_filters', 'single_photon_triggers',
					'exclude_wjets_greater_400', 'exclude_wjets_less_400',
					'one_ak15',
					'leading_fj250',
					'isoneG',
					'noextrab',
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
			sname = 'nominal' if systematic is None else systematic
			if systematic in weights.variations:
				weight = weights.weight(modifier=systematic)[cut]
			else:
				weight = weights.weight()[cut]
			#output['template'].fill(
			#	  region=region,
			#	  systematic=sname,
			#	  #recoil=normalize(u[region].r, cut),
			#	  #fjmass=normalize(leading_fj.pt, cut),
			#	  #TvsQCD=normalize(leading_fj.TvsQCD, cut),
			#	  weight=weight
			#)
			if systematic is None:
				variables = {
					'met':					met.pt,
					'metphi':				 met.phi,
					'ut' :				   u[region].r,
					'uphi' :				 u[region].phi,
					'mindphirecoil':		  ak.min(abs(u['sr'].delta_phi(j_clean.T)), axis=1,mask_identity=False),
					'minDphirecoil':		  abs(u[region].delta_phi(leading_fj.T)),
					'mupt':					 leading_mu.pt,
					'muphi':				 leading_mu.phi,
					'mueta':				 leading_mu.eta,
					'elept':					 leading_e.pt,
					'elephi':				 leading_e.phi,
					'eleeta':				 leading_e.eta,
					'phopt':					 leading_pho.pt,
					'phophi':				 leading_pho.phi,
					'phoeta':				 leading_pho.eta,
					'fjpt':					 leading_fj.pt,
					'fjphi':				 leading_fj.phi,
					'fjeta':				 leading_fj.eta,
					'jpt':					 leading_j.pt,
					'jphi':				 leading_j.phi,
					'jeta':				 leading_j.eta,
					#'jbtagLpt':					 j_btagvL.pt,
					#'jbtagLphi':				 j_btagvL.phi,
					#'jbtagLeta':				 j_btagvL.eta,
					'njbtagL':				 j_nbtagvL,
					'njets':				 j_nclean,
				}
				if region in mT:
					variables['mT']		   = mT[region]
				if 'z' in region:
					if 'e' in region:
						variables['dielept']	 = diele.pt
						variables['dielephi']	 = diele.phi
						variables['dieleeta']	 = diele.eta
					if 'z' in region:
						variables['dimupt']		= dimu.pt
						variables['dimuphi']	 = dimu.phi
						variables['dimueta']	 = dimu.eta

				for variable in output:
					if variable not in variables:
						continue
					normalized_variable = {variable: normalize(variables[variable],cut)}
					output[variable].fill(
						region=region,
						**normalized_variable,
						weight=weight,
					)
				output['TvsQCD'].fill(
					  region=region,
					  TvsQCD=normalize(leading_fj.TvsQCD, cut),
					  weight=weight
				)
				output['probT'].fill(
					  region=region,
					  probT=normalize(leading_fj.probT, cut),
					  weight=weight
				)
				output['probQCD'].fill(
					  region=region,
					  probQCD=normalize(leading_fj.probQCD, cut),
					  weight=weight
				)


		shift_name = None
		if shift_name is None:
			systematics = [None] + list(weights.variations)
		else:
			systematics = [shift_name]
			
		for region in regions:
			if region not in selected_regions: continue


			###
			# Adding recoil and minDPhi requirements
			###

			selection.add('recoil_'+region, (u[region].r>350))
			selection.add('mindphi_'+region, (ak.min(abs(u['sr'].delta_phi(j_clean.T)), axis=1, mask_identity=False) > 0.5))
			#selection.add('minDphi_'+region, (ak.min(abs(u[region].delta_phi(leading_fj.T)), axis=1, mask_identity=False) > 1.5))
			selection.add('minDphi_'+region, (abs(u[region].delta_phi(leading_fj.T)) > 1.5))
			regions[region].insert(6, 'recoil_'+region)
			regions[region].insert(7, 'mindphi_'+region)
			regions[region].insert(8, 'minDphi_'+region)

			for systematic in systematics:
				if isData and systematic is not None:
					continue
				fill(region, systematic)

			#########
			## Store cutflow
			#########
			vcut=np.zeros(len(events), dtype=np.int)
			output['cutflow'].fill(region=region,cutname='Initial', cutflow=vcut, weight=weights.weight())
			cuts = regions[region]
			allcuts = set()
			for i, icut in enumerate(cuts):
				allcuts.add(icut)
				jcut = selection.all(*allcuts)
				vcut = (i+1)*jcut
				output['cutflow'].fill(region=region,cutname=icut, cutflow=vcut, weight=weights.weight()*jcut)

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
