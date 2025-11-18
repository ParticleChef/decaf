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
			'nocr' :('WtoLNu-2Jets','DYto2L', 'Muon'),
			'pucr' :('WtoLNu-2Jets','DYto2L', 'Muon'),
			'allcr' :('WtoLNu-2Jets','DYto2L', 'Muon'),
		}
		

		self._met_triggers = { ## NOT UPDATED for Run3
			'2022pre': [
				'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
				'PFMETNoMu120_PFMHTNoMu120_IDTight'
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
			#'nfjets': hist.Hist(
			#	hist.axis.StrCategory([], name='region', growth=True),
			#	hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='nfjets', label='Number of AK15 Jets'),
			#	storage=hist.storage.Weight(),
			#),
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
            'j1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(38,250,1200, name='j1pt', label='AK4 Leading Jet $p_{T}$'),
                storage=hist.storage.Weight(),
            ),
            'j1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1eta', label='AK4 Leading Jet Eta'),
                storage=hist.storage.Weight(),
            ),
            'j1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1phi', label='AK4 Leading Jet Phi'),
                storage=hist.storage.Weight(),
            ),
            'fj1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(38,250,1200, name='fj1pt', label='AK15 Leading Jet Pt'),
                storage=hist.storage.Weight(),
            ),
            'fj1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1eta', label='AK15 Leading Jet Eta'),
                storage=hist.storage.Weight(),
            ),
            'fj1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1phi', label='AK15 Leading Jet Phi'),
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
		#isLooseElectron = self._ids['isLooseElectron'] 
		#isTightElectron = self._ids['isTightElectron'] 
		#isLoosePhoton   = self._ids['isLoosePhoton']   
		#isTightPhoton   = self._ids['isTightPhoton']   
		#isGoodAK4	   = self._ids['isGoodAK4']	   
		#isGoodAK15	= self._ids['isGoodAK15']	
		
		#deepflavWPs = self._common['btagWPs']['deepflav'][self._year]
		#deepcsvWPs = self._common['btagWPs']['deepcsv'][self._year]


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
		mu['id_sf'] = ak.where(
			mu.istight, 
			get_mu_tight_id_sf(self._year, abs(mu.eta), mu.pt), 
			mu.id_sf
		)
		mu['iso_sf'] = ak.where(
			mu.istight, 
			get_mu_tight_iso_sf(self._year, abs(mu.eta), mu.pt), 
			ak.ones_like(mu.pt) #mu.iso_sf
		)
		mu['highpt_id'] = ak.where(
			mu.istight, 
			get_mu_highpt_id_sf(self._year, abs(mu.eta), mu.pt), 
			ak.ones_like(mu.pt)
		)
		mu['hlt_id'] = ak.where(
			mu.istight, 
			get_mu_hlt_sf(self._year, abs(mu.eta), mu.pt), 
			ak.ones_like(mu.pt)
		)
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
		print('mu tight: ', (mu_ntight==1))
		print('mu loose: ', (mu_nloose==0))
		print('and: ', (mu_ntight==1)&(mu_nloose==0))
		
		# define leading mu
		leading_mu = ak.firsts(mu_tight)
		
		# define second mu for Z->mumu
		second_mu = ak.pad_none(mu_tight, target=2)[:,1]
		
		#dimu = ak.cartesian({"mu1": leading_mu, "mu2": second_mu}, nested=True)
		dimu = leading_mu + second_mu
		
		#dimu_mass = (dimu.mu1 + dimu.mu2).mass
		dimu_mass = dimu.mass


		j = events.Jet
		#j['isgood'] = isGoodAK4(j, self._year)
		#j['isclean'] = (
		#	ak.all(j.metric_table(mu_loose) > 0.4, axis=2)
		#	& ak.all(j.metric_table(e_loose) > 0.4, axis=2)
		#)
		#j['isdflvL'] = (j.btagDeepFlavB>deepflavWPs['loose'])
		#j['T'] = ak.zip(
		#	{
		#		"r": j.pt,
		#		"phi": j.phi,
		#	},
		#	with_name="PolarTwoVector",
		#	behavior=vector.behavior,
		#)
		#j_good = j[j.isgood]
		#j_clean = j_good[j_good.isclean]
		#j_dflvL = j_iso[j_clean.isdflvL]
		#leading_j = ak.firsts(j_clean)

		#j_ntot=ak.num(j, axis=1)
		#j_ngood=ak.num(j_good, axis=1)
		#j_nclean=ak.num(j_clean, axis=1)
		#j_ndflvL=ak.num(j_dflvL, axis=1)

		###
		# Calculate recoil and transverse mass
		###

		u = {
			'nocr'  : met+leading_mu.T,
			'pucr'  : met+leading_mu.T+second_mu.T,
			'allcr'  : met+leading_mu.T+second_mu.T,
		}

		mT = {
			'nocr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
			'pucr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
			'allcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
		}

		###
		#Calculating weights
		###
		if not isData:
			
			#gen = events.GenPart

			#gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			#gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			#gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			#genTops = gen[gen.isTop]
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
				'nocr' :  np.ones(len(events), dtype='float'),
				'pucr' : get_pu_weight(self._year, events.Pileup.nTrueInt),
				'allcr': get_pu_weight(self._year, events.Pileup.nTrueInt),
			}

			###
			# Trigger efficiency weight
			###
			#trig = {
			#	'nocr':   get_met_trig_weight(self._year, met.pt),
			#	'pucr': get_met_trig_weight(self._year, u['wmcr'].r),
			#}

			### 
			# Calculating electron and muon ID weights
			###
			ids ={
				'nocr':  np.ones(len(events), dtype='float'),
				'pucr':  np.ones(len(events), dtype='float'),
				'allcr': leading_mu.highpt_id#get_mu_highpt_id_sf(self._year, abs(leading_mu.eta), leading_mu.pt),
			}
		   
			###
			# Reconstruction weights for electrons
			###									   
			reco = {
				'nocr': np.ones(len(events), dtype='float'),
				'pucr': np.ones(len(events), dtype='float'),
				'allcr': np.ones(len(events), dtype='float'),
			}

			###
			# Isolation weights for muons
			###
			isolation = {
				'nocr':  np.ones(len(events), dtype='float'),
				'pucr':  np.ones(len(events), dtype='float'),
				'allcr': leading_mu.iso_sf#get_mu_tight_iso_sf(self._year, abs(leading_mu.eta), leading_mu.pt),
			}

			###
			# HLT weights for muons
			###
			hlt= {
				'nocr':  np.ones(len(events), dtype='float'),
				'pucr':  np.ones(len(events), dtype='float'),
				'allcr': leading_mu.hlt_id#get_mu_hlt_sf(self._year, abs(leading_mu.eta), leading_mu.pt),
			}

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
			weights.add('pileup',pu[region])
			weights.add('ids', ids[region])
			weights.add('isolation', isolation[region])
			weights.add('hlt', hlt[region])
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
		#selection.add('met_triggers', triggers)

		triggers = np.zeros(len(events), dtype='bool')
		for path in self._singlemuon_triggers[self._year]:
			if path not in events.HLT.fields:
				continue
			triggers = triggers | events.HLT[path]
		selection.add('single_muon_triggers', ak.to_numpy(triggers))

		

		selection.add('isoneM', ((mu_ntight==1)&(mu_nloose==0)))

		selection.add('tiMuonpt', (leading_mu.pt>53))
		
		#selection.add('one_ak4', (j_nclean>0))
		#selection.add('noextrab', (j_ndflvL==0))
		#selection.add('extrab', (j_ndflvL>0))


		#selection.add('met120',(met.pt<120))
		#selection.add('met150',(met.pt>150))

		regions = {
			'nocr': [
					'lumimask',
					'met_filters',
					'single_muon_triggers',
					'isoneM',
					'tiMuonpt'
			],
			'pucr': [
					'lumimask',
					'met_filters',
					'single_muon_triggers',
					'isoneM',
					'tiMuonpt'
			],
			'allcr': [
					'lumimask',
					'met_filters',
					'single_muon_triggers',
					'isoneM',
					'tiMuonpt'
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
			output['template'].fill(
				  region=region,
				  systematic=sname,
				  #recoil=normalize(u[region].r, cut),
				  #fjmass=normalize(leading_fj.pt, cut),
				  #TvsQCD=normalize(leading_fj.TvsQCD, cut),
				  weight=weight
			)
			if systematic is None:
				variables = {
					'met':					met.pt,
					'metphi':				 met.phi,
					'mupt':				  leading_mu.pt,
					'muphi':				 leading_mu.phi,
					'mueta':				 leading_mu.eta,
				}
				if region in mT:
					variables['mT']		   = mT[region]
				#if 'e' in region:
				#	if 'z' in region:
				#		variables['l1pt']	  = diele.pt
				#		variables['l1phi']	 = diele.phi
				#		variables['l1eta']	 = diele.eta
				#	else:
				#		variables['l1pt']	  = leading_e.pt
				#		variables['l1phi']	 = leading_e.phi
				#		variables['l1eta']	 = leading_e.eta
				if 'm' in region:
					if 'z' in region:
						variables['l1pt']	  = dimu.pt
						variables['l1phi']	 = dimu.phi
						variables['l1eta']	 = dimu.eta
					else:
						variables['l1pt']	  = leading_mu.pt
						variables['l1phi']	 = leading_mu.phi
						variables['l1eta']	 = leading_mu.eta

				for variable in output:
					if variable not in variables:
						continue
					normalized_variable = {variable: normalize(variables[variable],cut)}
					output[variable].fill(
						region=region,
						**normalized_variable,
						weight=weight,
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
				print(i, ", ", icut, ": ", events.run[jcut], " len: ", len(events.run[jcut]))

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
