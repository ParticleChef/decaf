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
import correctionlib

def update(events, collections):
	"""Return a shallow copy of events array with some collections swapped out"""
	out = events
	for name, value in collections.items():
		out = ak.with_field(out, value, name)
	return out

###
# Helpers for weight up/down variations (top-tag SF uncertainty study,
# toptag/docs/toptag_sf_uncertainty_plan.txt B2). Evaluated here instead of in
# utils/corrections.py so the shared corrections.coffea stays unchanged.
###
_PU_KEYS = {
	'2022pre' : 'Collisions2022_355100_357900_eraBCD_GoldenJson',
	'2022post': 'Collisions2022_359022_362760_eraEFG_GoldenJson',
	'2023pre' : 'Collisions2023_366403_369802_eraBC_GoldenJson',
	'2023post': 'Collisions2023_369803_370790_eraD_GoldenJson',
}
_CSET_CACHE = {}
def _cset(path):
	if path not in _CSET_CACHE:
		_CSET_CACHE[path] = correctionlib.CorrectionSet.from_file(path)
	return _CSET_CACHE[path]

def finite_or(x, fallback):
	"""numpy float64 copy of x with non-finite entries replaced by fallback."""
	x = np.asarray(x, dtype=np.float64)
	return np.where(np.isfinite(x), x, fallback)

def get_pu_weight_var(year, trueint, var):
	"""Pileup weight, var in {'nominal','up','down'} (absolute)."""
	trueint = np.asarray(ak.to_numpy(trueint), dtype=np.float64)
	return _cset('data/PUweight/'+year+'/puWeights.json.gz')[_PU_KEYS[year]].evaluate(trueint, var)

def get_mu_sf_var(year, key, eta, pt, var):
	"""Muon POG SF (muon_Z.json.gz) for jagged eta/pt, var in
	{'nominal','systup','systdown'}; same clamping as utils/corrections.py."""
	eta = ak.where(eta > 2.399, ak.full_like(eta, 2.399), eta)
	pt  = ak.where(pt < 15., ak.full_like(pt, 15.), pt)
	counts = ak.num(eta)
	flateta = np.asarray(ak.to_numpy(ak.flatten(eta)), dtype=np.float64)
	flatpt  = np.asarray(ak.to_numpy(ak.flatten(pt)), dtype=np.float64)
	w = _cset('data/MuonSF/'+year+'/muon_Z.json.gz')[key].evaluate(flateta, flatpt, var)
	return ak.unflatten(w, counts)

def nano_weight_table(events, name, width):
	"""(n_entries, table) for a jagged nanoAOD weight vector (PSWeight,
	LHEScaleWeight, LHEPdfWeight); table is [nevents, width] with NaN padding.
	None if the branch is absent (e.g. data or a sample without it)."""
	if name not in events.fields:
		return None
	arr = events[name]
	n = np.asarray(ak.to_numpy(ak.num(arr, axis=1)))
	tab = ak.to_numpy(ak.fill_none(ak.pad_none(arr, width, axis=1, clip=True), np.nan))
	return n, np.asarray(tab, dtype=np.float64)

def get_ps_weights(events):
	"""ISR and FSR (up, down) relative weights from PSWeight (4 entries:
	[0] ISR=2 FSR=1, [1] ISR=1 FSR=2, [2] ISR=0.5 FSR=1, [3] ISR=1 FSR=0.5).
	Events without 4 entries get 1."""
	nev = len(events)
	one = np.ones(nev)
	res = nano_weight_table(events, 'PSWeight', 4)
	if res is None:
		return (one, one), (one, one)
	n, t = res
	ok = (n == 4)
	col = lambda i: np.where(ok, finite_or(t[:, i], 1.), 1.)
	return (col(0), col(2)), (col(1), col(3))

def get_scale_weights(events):
	"""7-point muR/muF envelope (up=max, down=min, relative) from
	LHEScaleWeight; drops the anti-correlated (0.5,2)/(2,0.5) points.
	9 entries: [4] nominal, use 0,1,3,5,7,8. 8 entries (nominal omitted):
	use 0,1,3,4,6,7. Otherwise 1."""
	nev = len(events)
	up, down = np.ones(nev), np.ones(nev)
	res = nano_weight_table(events, 'LHEScaleWeight', 9)
	if res is None:
		return up, down
	n, t = res
	for nent, idx in ((9, [0, 1, 3, 5, 7, 8]), (8, [0, 1, 3, 4, 6, 7])):
		sel = (n == nent)
		if not np.any(sel):
			continue
		v = t[sel][:, idx]
		v = np.where(np.isfinite(v), v, 1.)
		up[sel] = np.max(v, axis=1)
		down[sel] = np.min(v, axis=1)
	return up, down

def get_pdf_weights(events):
	"""PDF and alpha_S (up, down) relative weights from LHEPdfWeight, Hessian
	sets (CMS Run 3 default NNPDF3.1 325300: [0] central, [1..100]
	eigenvectors, [101],[102] alpha_S). delta = sqrt(sum_i (w_i-w_0)^2)/w_0,
	up/down = 1 +- delta. alpha_S: envelope of the two alpha_S members.
	101/103 -> 100 members, 31/33 -> 30 members; other lengths get 1."""
	nev = len(events)
	pdf_up, pdf_dn = np.ones(nev), np.ones(nev)
	as_up, as_dn = np.ones(nev), np.ones(nev)
	res = nano_weight_table(events, 'LHEPdfWeight', 103)
	if res is None:
		return (pdf_up, pdf_dn), (as_up, as_dn)
	n, t = res
	for nent, nmem, has_as in ((101, 100, False), (103, 100, True), (31, 30, False), (33, 30, True)):
		sel = (n == nent)
		if not np.any(sel):
			continue
		w0 = t[sel, 0]
		good = np.isfinite(w0) & (w0 > 0)
		w0s = np.where(good, w0, 1.)
		mem = t[sel, 1:nmem+1]
		d = np.where(np.isfinite(mem), mem - w0s[:, None], 0.)
		delta = np.where(good, np.sqrt(np.sum(d**2, axis=1)) / w0s, 0.)
		pdf_up[sel] = 1. + delta
		pdf_dn[sel] = np.maximum(1. - delta, 0.)
		if has_as:
			a = t[sel, nmem+1:nmem+3] / w0s[:, None]
			a = np.where(np.isfinite(a) & good[:, None], a, 1.)
			as_up[sel] = np.max(a, axis=1)
			as_dn[sel] = np.min(a, axis=1)
	return (pdf_up, pdf_dn), (as_up, as_dn)

class AnalysisProcessor(processor.ProcessorABC):

	lumis = { 
		#https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiRecommendationsRun3
		'2022pre' : 7.99,
		'2022post': 26.68,
		'2023pre' : 17.96,
		'2023post': 9.68,
	}

	lumiMasks = {
		'2022pre' : LumiMask("data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json"),
		'2022post': LumiMask("data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json"),
		'2023pre' : LumiMask("data/lumiMask/Cert_Collisions2023_366442_370790_Golden.json"),
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

		# Top-tag SF uncertainty study: only the regions the SF is measured in
		# (gcr -> SF_QCD, tmcr + tecr -> SF_Top) are filled; every other region
		# is skipped to save time and output size. See
		# toptag/docs/toptag_sf_uncertainty_plan.txt (B4).
		self._syst_regions = ('gcr', 'tmcr', 'tecr')

		self._TvsQCDwp = {
			'2022pre':  0.60,
			'2022post': 0.60,
			'2023pre':  0.53,
			'2023post': 0.53,
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
		self._photon_triggers = { ## Name is same 22/23
			'2022pre':
				[
				'Photon200'
			],
		}
		self._corrections = corrections
		self._ids = ids
		self._common = common

		self.make_output = lambda: {
			'sumw': 0.,
			# Only what the SF chain reads (toptag/studies_unc): fjpt_bin, fjpt_cat,
			# TvsQCD_cat (ROC) and sumw. Control plots and the cutflow come from
			# topstudy_run.py. n_qinfj/n_binfj: Integer(0,5); values >= 5 go to
			# overflow, which slice(4, None, sum) ("gt3") and sum still include.
			# fjpt_cat / fjpt_bin carry a 'systematic' axis ("nominal" + every
			# weights.variations entry; data: "nominal" only) for the top-tag
			# SF uncertainty study (toptag/docs/toptag_sf_uncertainty_plan.txt).
			'fjpt_cat': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.StrCategory([], name='systematic', growth=True),
				hist.axis.Variable([250, 325, 400, 600, 2000], name='fjpt', label='AK15 Leading Jet $p_{T}$', flow=True),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				hist.axis.Integer(0,5, name='n_qinfj', label='Number of q from top in AK15 Jet'),
				hist.axis.Integer(0,5, name='n_binfj', label='Number of b flavor in AK15 Jet'),
				storage=hist.storage.Weight(),
			),
			'fjpt_bin': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.StrCategory([], name='systematic', growth=True),
				#hist.axis.Regular(38,250,1200, name='fjpt', label='AK15 Leading Jet $p_{T}$'),
				hist.axis.Variable([250, 325, 400, 600, 2000], name='fjpt_bin', label='AK15 Leading Jet $p_{T}$', flow=True),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'TvsQCD_cat': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(15,0,1, name='TvsQCD', label='TvsQCD'),
				hist.axis.Integer(0,5, name='n_qinfj', label='Number of q from top in AK15 Jet'),
				hist.axis.Integer(0,5, name='n_binfj', label='Number of b flavor in AK15 Jet'),
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
		caller = dataset.split("____")[0]

		selected_regions = []
		for region, samples in self._samples.items():
			for sample in samples:
				if sample not in dataset: continue
				selected_regions.append(region)
		# keep only the SF measurement regions (gcr, tmcr, tecr)
		selected_regions = [r for r in selected_regions if r in self._syst_regions]

		isData = not hasattr(events, "genWeight")
		
		################
		## Initialize selection, weights and output histogram
		################
		selection = PackedSelection(dtype="uint64")
		weights = Weights(len(events), storeIndividual=True)
		output = self.make_output()

		if not isData:
			genw = ak.to_numpy(events.genWeight).astype(np.float64)
			weights.add("genw", genw)

			if shift_name is None:
				output["sumw"] = np.sum(genw, dtype=np.float64)
		#if shift_name is None and not isData:
		#	output['sumw'] = ak.sum(events.genWeight)

		###
		#Getting corrections, ids from .coffea files
		###

		get_pu_weight			= self._corrections['get_pu_weight']	
		get_mu_loose_id_sf	   = self._corrections['get_mu_loose_id_sf']
		get_mu_tight_id_sf	   = self._corrections['get_mu_tight_id_sf']
		get_mu_loose_iso_sf	  = self._corrections['get_mu_loose_iso_sf']
		get_mu_tight_iso_sf	  = self._corrections['get_mu_tight_iso_sf']
		get_met_xy_correction	= self._corrections['get_met_xy_correction']
		get_ele_loose_id_sf	  = self._corrections['get_ele_loose_id_sf']
		get_ele_tight_id_sf	  = self._corrections['get_ele_tight_id_sf']
		get_ele_reco_sf_below20	  = self._corrections['get_ele_reco_sf_below20']
		get_ele_reco_sf_20to75	  = self._corrections['get_ele_reco_sf_20to75']
		get_ele_reco_sf_Above75	  = self._corrections['get_ele_reco_sf_Above75']
		get_photon_id_sf	  = self._corrections['get_photon_id_sf']
		get_met_trig_weight	  = self._corrections['get_met_trig_weight']
		get_ele_trig_weight	  = self._corrections['get_ele_trig_weight']
		get_pho_trig_weight	  = self._corrections['get_pho_trig_weight']
		get_nlo_ewk_weight	   = self._corrections['get_nlo_ewk_weight']	
		#get_nnlo_nlo_weight	  = self._corrections['get_nnlo_nlo_weight'][self._year]
		get_btag_weight	  = self._corrections['get_btag_weight']
		get_ttbar_weight	 = self._corrections['get_ttbar_weight']
		#get_mu_rochester_sf	  = self._corrections['get_mu_rochester_sf'][self._year]

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
		mu['id_sf'] = ak.where(
			mu.isloose, 
			get_mu_loose_id_sf(self._year, abs(mu.eta), mu.pt), 
			ak.ones_like(mu.pt)
		)
		mu['iso_sf'] = ak.where(
			mu.isloose, 
			get_mu_loose_iso_sf(self._year, abs(mu.eta), mu.pt), 
			ak.ones_like(mu.pt)
		)
		mu['istight'] = isTightMuon(mu,self._year)
		mu['id_sf'] = ak.where(
			mu.istight, 
			get_mu_tight_id_sf(self._year, abs(mu.eta), mu.pt), 
			mu.id_sf
		)
		mu['iso_sf'] = ak.where(
			mu.istight, 
			get_mu_tight_iso_sf(self._year, abs(mu.eta), mu.pt), 
			mu.iso_sf #ak.ones_like(mu.pt) #mu.iso_sf
		)
		# up/down of the muon ID and ISO SFs (absolute, POG 'systup'/'systdown')
		if not isData:
			for sfname, lkey, tkey in (('id_sf', 'NUM_LooseID_DEN_TrackerMuons', 'NUM_TightID_DEN_TrackerMuons'),
			                           ('iso_sf', 'NUM_LoosePFIso_DEN_LooseID', 'NUM_TightPFIso_DEN_TightID')):
				for sfvar, pogvar in (('up', 'systup'), ('down', 'systdown')):
					v = ak.where(mu.isloose, get_mu_sf_var(self._year, lkey, abs(mu.eta), mu.pt, pogvar), ak.ones_like(mu.pt))
					v = ak.where(mu.istight, get_mu_sf_var(self._year, tkey, abs(mu.eta), mu.pt, pogvar), v)
					mu[sfname+'_'+sfvar] = v
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
		sf_below20, sf_below20_up, sf_below20_down = get_ele_reco_sf_below20(self._year, abs(e.eta), e.pt, e.phi)
		sf_20to75,  sf_20to75_up,  sf_20to75_down  = get_ele_reco_sf_20to75( self._year, abs(e.eta), e.pt, e.phi)
		sf_Above75, sf_Above75_up, sf_Above75_down = get_ele_reco_sf_Above75(self._year, abs(e.eta), e.pt, e.phi)
		e['reco_sf'] = ak.where(
			e.pt < 20, 
			sf_below20, 
			ak.where(e.pt < 75, sf_20to75, sf_Above75)
		)
		sf_loose, sf_loose_up, sf_loose_down = get_ele_loose_id_sf(self._year, abs(e.eta), e.pt, e.phi)
		sf_tight, sf_tight_up, sf_tight_down = get_ele_tight_id_sf(self._year, abs(e.eta), e.pt, e.phi)
		eleHlt_loose, eleHlt_loose_up, eleHlt_loose_down = get_ele_trig_weight(self._year, abs(e.eta), e.pt, "Loose")
		eleHlt_tight, eleHlt_tight_up, eleHlt_tight_down = get_ele_trig_weight(self._year, abs(e.eta), e.pt, "Tight")
		e['isloose'] = isLooseElectron(e,self._year)
		e['id_sf'] = ak.where(
			e.isloose, 
			sf_loose, 
			ak.ones_like(e.pt)
		)
		e['hlt_sf'] = ak.where(
			e.isloose, 
			eleHlt_loose, 
			ak.ones_like(e.pt)
		)
		e['istight'] = isTightElectron(e,self._year)
		e['id_sf'] = ak.where(
			e.istight, 
			sf_tight, 
			e.id_sf #ak.ones_like(e.pt)
		)
		e['hlt_sf'] = ak.where(
			e.istight, 
			eleHlt_tight, 
			e.hlt_sf
		)
		# up/down of the electron RECO, ID and HLT SFs (absolute), same
		# loose -> tight override as the nominal
		for v in ('up', 'down'):
			below20, s20to75, above75 = {'up': (sf_below20_up, sf_20to75_up, sf_Above75_up),
			                             'down': (sf_below20_down, sf_20to75_down, sf_Above75_down)}[v]
			e['reco_sf_'+v] = ak.where(e.pt < 20, below20, ak.where(e.pt < 75, s20to75, above75))
			id_l, id_t = {'up': (sf_loose_up, sf_tight_up), 'down': (sf_loose_down, sf_tight_down)}[v]
			e['id_sf_'+v] = ak.where(e.istight, id_t, ak.where(e.isloose, id_l, ak.ones_like(e.pt)))
			hlt_l, hlt_t = {'up': (eleHlt_loose_up, eleHlt_tight_up), 'down': (eleHlt_loose_down, eleHlt_tight_down)}[v]
			e['hlt_sf_'+v] = ak.where(e.istight, hlt_t, ak.where(e.isloose, hlt_l, ak.ones_like(e.pt)))
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
		pho_sf_loose, pho_sf_loose_up, pho_sf_loose_down = get_photon_id_sf(self._year, 'Loose', abs(pho.eta), pho.pt, pho.phi)
		pho_sf_tight, pho_sf_tight_up, pho_sf_tight_down = get_photon_id_sf(self._year, 'Tight', abs(pho.eta), pho.pt, pho.phi)
		pho['isloose'] = isLoosePhoton(pho,self._year)
		pho['id_sf'] = ak.where(
			pho.isloose, 
			pho_sf_loose, 
			ak.ones_like(pho.pt)
		)
		pho['istight'] = isTightPhoton(pho,self._year)
		pho['id_sf'] = ak.where(
			pho.istight, 
			pho_sf_tight, 
			pho.id_sf #ak.ones_like(pho.pt)
		)
		pho['id_sf_up'] = ak.where(pho.istight, pho_sf_tight_up, ak.where(pho.isloose, pho_sf_loose_up, ak.ones_like(pho.pt)))
		pho['id_sf_down'] = ak.where(pho.istight, pho_sf_tight_down, ak.where(pho.isloose, pho_sf_loose_down, ak.ones_like(pho.pt)))
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
		rho_density = events.Rho.fixedGridRhoFastjetAll
		fjec_corr = get_fjec_correction(self._year, fj.pt, fj.eta, fj.phi, rho_density, fj.area, run, isData)
		fj['pt'] = fj.pt * fjec_corr
		fj['mass'] = fj.mass * fjec_corr
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
			'sr'	: u['sr'].delta_phi(leading_fj.T),
			'wmcr'  : u['wmcr'].delta_phi(leading_fj.T),
			'wecr'  : u['wecr'].delta_phi(leading_fj.T),
			'tmcr'  : u['tmcr'].delta_phi(leading_fj.T),
			'tecr'  : u['tecr'].delta_phi(leading_fj.T),
			'zmcr'  : u['zmcr'].delta_phi(leading_fj.T),
			'zecr'  : u['zecr'].delta_phi(leading_fj.T),
			'gcr'   : u['gcr'].delta_phi(leading_fj.T)
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
			#j['isbflav'] = j.hadronFlavour == 5
			#j_bflav = j[j.isbflav]
			#binfj = leading_fj.vec.delta_r(j_bflav)
			#n_binfj = ak.sum(binfj < 1.5, axis=1)
			n_binfj = ak.fill_none(leading_fj.nBHadrons, 0)
			
			gen = events.GenPart

			yourmotheris = gen.distinctParent.distinctParent.pdgId
			gen['isq'] = (abs(gen.pdgId) >= 1) & (abs(gen.pdgId) <= 5) & gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			p1 = ak.fill_none(gen.distinctParent.pdgId, 0)
			p2 = ak.fill_none(ak.fill_none(gen.distinctParent.distinctParent.pdgId, 0), 0)
			gen['isqfromt'] = gen.isq & (abs(p1) == 24) & (abs(p2) == 6)
			gen['isbfromt'] = gen.isq & (abs(gen.pdgId) == 5) & (abs(p1) == 6)
			q_from_t = gen[gen.isqfromt]
			b_from_t = gen[gen.isbfromt]

			b_top = b_from_t.distinctParentIdxG
			q_top = q_from_t.distinctParent.distinctParentIdxG
			b_from_t = b_from_t[ak.sum(b_top[:, :, None] == q_top[:, None, :], axis=2) >= 2]
			q_from_t = ak.concatenate([q_from_t, b_from_t], axis=1)

			qinfj = leading_fj.vec.delta_r(q_from_t)
			n_qinfj = ak.sum(qinfj <= 1.5, axis=1)
			#fj['qisinfj'] = ak.all(fj.vec.metric_table(q_from_t)<1.5, axis=2)
			#q_from_t_fj = fj[fj.qisinfj]
			#nq_from_t_fj = ak.num(q_from_t_fj, axis=1)

			gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			genTops = gen[gen.isTop]
			genW = gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))]
			gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isA'] = (abs(gen.pdgId)==22) & gen.hasFlags(['isPrompt', 'fromHardProcess', 'isLastCopy']) & (gen.status == 1)
			#
			genWs = gen[gen.isW] 
			genZs = gen[gen.isZ]
			genDYs = gen[gen.isZ&(gen.mass>30)]
			genAs = gen[gen.isA & (gen.pt > 100)]
			#
			nlo = np.ones(len(events), dtype='float')
			if('TTto' in dataset): 
				nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt) * get_ttbar_weight(genTops[:,1].pt))
			nlo_ewk = np.ones(len(events), dtype='float')
			if('WtoLNu' in dataset): 
				nlo_ewk = get_nlo_ewk_weight['w'](ak.max(genWs.pt, axis=1))
			elif('DY' in dataset): 
				nlo_ewk = get_nlo_ewk_weight['dy'](ak.max(genDYs.pt, axis=1))
			elif('Zto2Nu' in dataset): 
				nlo_ewk = get_nlo_ewk_weight['z'](ak.max(genZs.pt, axis=1))
			elif('GJ' in dataset):
				nlo_ewk = get_nlo_ewk_weight['a'](ak.max(genAs.pt, axis=1))

			GenIsoPho = events.GenIsolatedPhoton
			nGenIsoPho = ak.num(GenIsoPho, axis=1)
			#
			###
			# Calculate PU weight and systematic variations
			###
			pu = get_pu_weight(self._year, events.Pileup.nTrueInt)
			pu_up = get_pu_weight_var(self._year, events.Pileup.nTrueInt, 'up')
			pu_down = get_pu_weight_var(self._year, events.Pileup.nTrueInt, 'down')

			###
			# AK4 b-tagging weights
			###
			btagSF, \
			btagSFbc_correlatedUp, \
			btagSFbc_correlatedDown, \
			btagSFbc_uncorrelatedUp, \
			btagSFbc_uncorrelatedDown, \
			btagSFlight_correlatedUp, \
			btagSFlight_correlatedDown, \
			btagSFlight_uncorrelatedUp, \
			btagSFlight_uncorrelatedDown  = get_btag_weight('PNetUParT',self._year,'loose',caller).btag_weight(
				j_iso.pt,
				j_iso.eta,
				j_iso.hadronFlavour,
				j_iso.isbtagvL
			)


		###
		# Selections
		###

		lumimask = np.ones(len(events), dtype='bool')
		if isData:
			lumimask = AnalysisProcessor.lumiMasks[self._year](events.run, events.luminosityBlock)
		selection.add('lumimask', lumimask)

		met_filters =  np.ones(len(events), dtype='bool')
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
		if not isData:
			# Weights.add(name, w, wUp, wDown): wUp/wDown are absolute (coffea
			# divides by w where w != 0); variations appear as <name>Up/Down.
			pu_nom = finite_or(pu, 1.)
			weights.add('pileup', pu_nom, finite_or(pu_up, pu_nom), finite_or(pu_down, pu_nom))
			weights.add('nlo_ewk',nlo_ewk)
			# top-pT reweighting: up = applied twice, down = not applied
			# (nlo = 1 for non-ttbar samples, so no effect there)
			nlo_nom = finite_or(ak.to_numpy(nlo) if not isinstance(nlo, np.ndarray) else nlo, 1.)
			weights.add('topptreweighting', nlo_nom, nlo_nom**2, np.ones(len(events)))
			# parton shower, muR/muF scale envelope, PDF and alpha_S (relative
			# weights; 1 where the sample has no such branch)
			(isr_up, isr_down), (fsr_up, fsr_down) = get_ps_weights(events)
			one = np.ones(len(events))
			weights.add('ISR', one, isr_up, isr_down)
			weights.add('FSR', one, fsr_up, fsr_down)
			scale_up, scale_down = get_scale_weights(events)
			weights.add('QCDscale', one, scale_up, scale_down)
			(pdf_up, pdf_down), (as_up, as_down) = get_pdf_weights(events)
			weights.add('PDF', one, pdf_up, pdf_down)
			weights.add('alphaS', one, as_up, as_down)
			weights.add('btagSF',btagSF)
			# b-tag variations as ratios to the nominal btagSF. btag_weight()
			# returns btagSF = 0 when an untagged jet has SF*eff >= 1 (data
			# efficiency clipped to 1); a plain Up/btagSF is then 0/0 = NaN and
			# poisons whole histogram bins of every variation. Those events have
			# zero nominal weight, so their variation ratio is set to 1.
			btag_nom = np.asarray(ak.to_numpy(btagSF), dtype=np.float64)
			btag_ok = np.isfinite(btag_nom) & (btag_nom > 0)
			btag_den = np.where(btag_ok, btag_nom, 1.)
			def btag_ratio(var):
				r = np.asarray(ak.to_numpy(var), dtype=np.float64) / btag_den
				return np.where(btag_ok & np.isfinite(r), r, 1.)
			weights.add('btagSFbc_correlated',np.ones(len(events), dtype='float'), btag_ratio(btagSFbc_correlatedUp), btag_ratio(btagSFbc_correlatedDown))
			weights.add('btagSFbc_uncorrelated',np.ones(len(events), dtype='float'), btag_ratio(btagSFbc_uncorrelatedUp), btag_ratio(btagSFbc_uncorrelatedDown))
			weights.add('btagSFlight_correlated',np.ones(len(events), dtype='float'), btag_ratio(btagSFlight_correlatedUp), btag_ratio(btagSFlight_correlatedDown))
			weights.add('btagSFlight_uncorrelated',np.ones(len(events), dtype='float'), btag_ratio(btagSFlight_uncorrelatedUp), btag_ratio(btagSFlight_uncorrelatedDown))

			### 
			# Region masks 
			###
			mask_iszeroL = ak.to_numpy((e_nloose==0)&(mu_nloose==0)&(pho_nloose==0))
			mask_isoneE = ak.to_numpy((e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(pho_nloose==0))
			mask_isoneM = ak.to_numpy((e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(pho_nloose==0))
			mask_istwoE = ak.to_numpy((e_nloose==2)&(mu_nloose==0)&(pho_nloose==0))
			mask_istwoM = ak.to_numpy((e_nloose==0)&(mu_nloose==2)&(pho_nloose==0))
			mask_isoneG = ak.to_numpy((e_nloose==0)&(mu_nloose==0)&(pho_nloose==1)&(pho_ntight==1))

			### 
			# Object SFs, one weight per object type so that electron, muon and
			# photon SFs vary independently. Each event gets the product over
			# its region's objects (masks are disjoint); the nominal product of
			# these weights equals the former 'ids'*'reco'*'iso'*'trig'.
			###
			nev = len(events)
			def obj_prod(objs, field):
				return finite_or(ak.to_numpy(ak.fill_none(ak.prod(objs[field], axis=1), 1.0)), 1.)
			def add_obj_sf(name, pieces):
				"""pieces: list of (mask, objects, field). Builds nominal/up/down
				from field, field+'_up', field+'_down' on the masked events."""
				w = {s: np.ones(nev) for s in ('', '_up', '_down')}
				for mask, objs, field in pieces:
					for s in w:
						w[s][mask] = obj_prod(objs, field+s)[mask]
				weights.add(name, w[''], w['_up'], w['_down'])

			# ID: electron (1e: tight, 2e: loose), muon (1mu: tight, 2mu: loose), photon (tight)
			add_obj_sf('ids_e',   [(mask_isoneE, e_tight, 'id_sf'), (mask_istwoE, e_loose, 'id_sf')])
			add_obj_sf('ids_mu',  [(mask_isoneM, mu_tight, 'id_sf'), (mask_istwoM, mu_loose, 'id_sf')])
			add_obj_sf('ids_pho', [(mask_isoneG, pho_tight, 'id_sf')])
			# RECO: electron
			add_obj_sf('reco_e',  [(mask_isoneE, e_tight, 'reco_sf'), (mask_istwoE, e_loose, 'reco_sf')])
			# ISO: muon
			add_obj_sf('iso_mu',  [(mask_isoneM, mu_tight, 'iso_sf'), (mask_istwoM, mu_loose, 'iso_sf')])
			# HLT: electron
			add_obj_sf('trig_e',  [(mask_isoneE, e_tight, 'hlt_sf'), (mask_istwoE, e_loose, 'hlt_sf')])

			###
			# Trigger efficiency weight: MET trigger (0l, 1mu, 2mu) and photon trigger (1g)
			###
			met_sf_0,  met_sf_0_up,  met_sf_0_down  = get_met_trig_weight(self._year, u['sr'].r)
			met_sf_m,  met_sf_m_up,  met_sf_m_down  = get_met_trig_weight(self._year, u['wmcr'].r)
			met_sf_mm, met_sf_mm_up, met_sf_mm_down = get_met_trig_weight(self._year, u['zmcr'].r)
			pho_sf, pho_sf_up, pho_sf_down = get_pho_trig_weight(self._year, ak.firsts(pho_tight).pt)
			def ev_arr(x):
				return finite_or(ak.to_numpy(ak.fill_none(x, 1.0)), 1.)

			trg_met = {s: np.ones(nev) for s in ('', '_up', '_down')}
			for mask, trio in ((mask_iszeroL, (met_sf_0,  met_sf_0_up,  met_sf_0_down)),
			                   (mask_isoneM,  (met_sf_m,  met_sf_m_up,  met_sf_m_down)),
			                   (mask_istwoM,  (met_sf_mm, met_sf_mm_up, met_sf_mm_down))):
				for s, x in zip(('', '_up', '_down'), trio):
					trg_met[s][mask] = ev_arr(x)[mask]
			weights.add('trig_met', trg_met[''], trg_met['_up'], trg_met['_down'])

			trg_pho = {s: np.ones(nev) for s in ('', '_up', '_down')}
			for s, x in zip(('', '_up', '_down'), (pho_sf, pho_sf_up, pho_sf_down)):
				trg_pho[s][mask_isoneG] = ev_arr(x)[mask_isoneG]
			weights.add('trig_pho', trg_pho[''], trg_pho['_up'], trg_pho['_down'])




		def normalize(val, cut,fillvalue=np.nan):
			if cut is None:
				ar = ak.to_numpy(ak.fill_none(val, fillvalue))
				return ar
			else:
				ar = ak.to_numpy(ak.fill_none(val[cut], fillvalue))
				return ar
			#if cut is None:
			#	ar = ak.to_numpy(ak.fill_none(val, np.nan))
			#	return ar
			#else:
			#	ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
			#	return ar
				
		def fill_region(region, systematics):
			# cut and the per-event arrays are the same for every systematic, so
			# they are computed once per region; only the weight changes.
			cut = selection.all(*regions[region])
			fjpt   = normalize(leading_fj.pt, cut)
			tvsqcd = normalize(leading_fj.TvsQCD, cut)
			if not isData:
				nq = normalize(n_qinfj, cut, fillvalue=-1)
				nb = normalize(n_binfj, cut, fillvalue=-1)
			for systematic in systematics:
				# data: "nominal" only
				if isData and systematic is not None:
					continue
				sname = 'nominal' if systematic is None else systematic
				if systematic in weights.variations:
					weight = weights.weight(modifier=systematic)[cut]
				else:
					weight = weights.weight()[cut]
				output['fjpt_bin'].fill(
					  region=region,
					  systematic=sname,
					  fjpt_bin=fjpt,
					  TvsQCD=tvsqcd,
					  weight=weight
				)
				if not isData:
					output['fjpt_cat'].fill(
						  region=region,
						  systematic=sname,
						  fjpt=fjpt,
						  TvsQCD=tvsqcd,
						  n_qinfj=nq,
						  n_binfj=nb,
						  weight=weight
					)
					if systematic is None:
						output['TvsQCD_cat'].fill(
							  region=region,
							  TvsQCD=tvsqcd,
							  n_qinfj=nq,
							  n_binfj=nb,
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

			fill_region(region, systematics)

		scale = 1
		if self._xsec[dataset]!= -1: 
			#scale = 1.0#self._lumi*self._xsec[dataset]
			scale = self._lumi*self._xsec[dataset]

		for key in output:
			if key=='sumw': 
				continue
			output[key] *= scale #* skimrate[dataset.split('____')]
				
		return output

	def postprocess(self, accumulator):

		return accumulator

if __name__ == '__main__':
	# `python3 processors/topstudy_unc.py` puts processors/ (not analysis/) on
	# sys.path; corrections_cached.coffea needs `utils.cset_cache` importable.
	import os, sys
	_ANALYSIS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	if _ANALYSIS_DIR not in sys.path:
		sys.path.insert(0, _ANALYSIS_DIR)

	parser = OptionParser()
	parser.add_option('-y', '--year', help='year', dest='year')
	parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
	parser.add_option('-n', '--name', help='name', dest='name')
	parser.add_option('-c', '--corrections', help='corrections file (default data/corrections.coffea)', dest='corrections', default='data/corrections.coffea')
	(options, args) = parser.parse_args()


	with gzip.open('metadata/'+options.metadata+'.json.gz') as fin:
		samplefiles = json.load(fin)
		xsec = {k: v['xs'] for k,v in samplefiles.items()}

	corrections = load(options.corrections)
	ids		 = load('data/ids.coffea')
	common	  = load('data/common.coffea')

	processor_instance=AnalysisProcessor(year=options.year,
										 xsec=xsec,
										 corrections=corrections,
										 ids=ids,
										 common=common)

	save(processor_instance, 'data/hadmonotop'+options.name+'.processor')
