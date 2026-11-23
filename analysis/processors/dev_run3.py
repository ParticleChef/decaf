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
# Helpers for weight up/down variations and JER smearing. Copied from
# processors/topstudy_unc.py (keep the two in sync). Evaluated here instead of
# in utils/corrections.py so the shared corrections.coffea stays unchanged.
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

###
# JER smearing (MC). Evaluated here like the helpers above; the resolution and
# scale factor come from the same jet_jerc / fatJet_jerc JSON files as the JEC
# (AK15 uses the AK8PFPuppi keys, as the JEC does).
###
_JER_TAG = {
	'2022pre' : 'Summer22_22Sep2023_JRV1',
	'2022post': 'Summer22EE_22Sep2023_JRV1',
	'2023pre' : 'Summer23Prompt23_RunCv1234_JRV1',	# full 2023 C era (also RunCv123, RunCv4 in the file)
	'2023post': 'Summer23BPixPrompt23_RunD_JRV1',
}
# Upper Rho edge of the PtResolution binning. Outside the binning (Rho or
# |eta| > 5.191) the JSON returns a resolution of 1.0 (100%), so the inputs
# are clamped into range before the lookup.
_JER_RHO_MAX = {'2022pre': 60.35, '2022post': 60.71, '2023pre': 52.2, '2023post': 52.05}
_JER_ETA_MAX = 5.191

def jer_factor_flat(year, pt, eta, rho, genpt, gendr, var, fat, seed):
	"""JER smearing factor for flat numpy arrays (one entry per jet), hybrid
	method; var in {'nom','up','down'} selects the JER scale factor s.
	  gen-matched (genpt > 0, dR < R/2, |pt-genpt| < 3*sigma*pt):
	      1 + (s-1)*(pt-genpt)/pt
	  otherwise: 1 + N(0,sigma)*sqrt(max(s^2-1, 0))
	sigma = relative pt resolution; pt = JEC-corrected pt. The factor is
	clipped at 0. The random numbers depend only on seed and the number of
	jets, so 'nom', 'up' and 'down' use the same ones."""
	fname, algo, cone = ('fatJet_jerc', 'AK8PFPuppi', 1.5) if fat else ('jet_jerc', 'AK4PFPuppi', 0.4)
	cset = _cset('data/JetMETCorr/'+year+'/'+fname+'.json.gz')
	pt, eta, rho, genpt, gendr = (np.asarray(x, dtype=np.float64) for x in (pt, eta, rho, genpt, gendr))
	eta_c = np.clip(eta, -_JER_ETA_MAX + 1e-3, _JER_ETA_MAX - 1e-3)
	rho_c = np.clip(rho, 0., _JER_RHO_MAX[year] - 1e-2)
	sigma = cset[_JER_TAG[year]+'_MC_PtResolution_'+algo].evaluate(eta_c, pt, rho_c)
	sf = cset[_JER_TAG[year]+'_MC_ScaleFactor_'+algo].evaluate(eta_c, pt, var)
	safe_pt = np.where(pt > 0, pt, 1.)
	matched = (genpt > 0) & (gendr < cone/2.) & (np.abs(pt - genpt) < 3.*sigma*pt)
	gaus = np.random.default_rng(seed).standard_normal(len(pt))
	scaling = 1. + (sf - 1.)*(pt - genpt)/safe_pt
	stochastic = 1. + gaus*sigma*np.sqrt(np.maximum(sf*sf - 1., 0.))
	factor = np.where(matched, scaling, stochastic)
	return np.where(np.isfinite(factor), np.maximum(factor, 0.), 1.)

def get_jer_factor(year, pt, eta, rho, genpt, gendr, var, fat, seed):
	"""jer_factor_flat for jagged per-jet arrays (rho is per event)."""
	counts = ak.num(pt)
	rho, _ = ak.broadcast_arrays(rho, pt)
	flat = lambda x: np.asarray(ak.to_numpy(ak.flatten(x)), dtype=np.float64)
	factor = jer_factor_flat(year, flat(pt), flat(eta), flat(rho), flat(genpt), flat(gendr), var, fat, seed)
	return ak.unflatten(factor, counts)

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
		# False: MC jets (AK15, AK4) get the nominal JER smearing and the
		# JERUp/JERDown passes run. True: no smearing at all, JES passes only.
		self._skipJER = False
		# Jet energy shifts run as extra passes over MC (see process()).
		# shift name -> (JES var for get_jec_correction / get_fjec_correction,
		#                JER var for get_jer_factor)
		self._jet_shifts = {'JESUp': ('up', 'nom'), 'JESDown': ('down', 'nom')}
		if not self._skipJER:
			self._jet_shifts.update({'JERUp': ('nominal', 'up'), 'JERDown': ('nominal', 'down')})

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
		self._ResolvedTvsQCDwp = { # ResolvedTvsQCD > 0.973943293094635
			'2022pre': 0.973943293094635,
			'2022post': 0.973943293094635,
			'2023pre': 0.973943293094635,
			'2023post': 0.973943293094635,
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
				# The only histogram with a 'systematic' axis: "nominal", every
				# weights.variations entry and JESUp/Down, JERUp/Down (data:
				# "nominal" only). Same recoil binning as 'ut'.
				hist.axis.Regular(175,250,2000, name='recoil', label='$U_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
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
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'minDphirecoil': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(30,0,3.5, name='minDphirecoil', label='Min |dPhi(Recoil, leading AK15s)|'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'ut': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(175,250,2000, name='ut', label='$U_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				#hist.axis.Variable([0, self._ResolvedTvsQCDwp[self._year], 1], name='ResolvedTvsQCD', label='ResolvedTvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'uphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='uphi', label='U_{phi}'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'mupt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='mupt', label='Leading Muon $p_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'mueta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='mueta', label='Leading Muon Eta'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'muphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='muphi', label='Leading Muon Phi'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'elept': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='elept', label='Leading Electron $p_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'eleeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='eleeta', label='Leading Electron Eta'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'elephi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='elephi', label='Leading Electron Phi'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'phopt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='phopt', label='Leading Photon $p_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'phoeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='phoeta', label='Leading Photon Eta'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'phophi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='phophi', label='Leading Photon Phi'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'dielept': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='dielept', label='DiElectron $p_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'dieleeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='dieleeta', label='DiElectron Eta'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'dielephi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='dielephi', label='DiElectron Phi'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'dimupt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='dimupt', label='DiMuon $p_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'dimueta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(48,-2.4,2.4, name='dimueta', label='DiMuon Eta'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'dimuphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(64,-3.2,3.2, name='dimuphi', label='DiMuon Phi'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'jpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,2000, name='jpt', label='AK4 Leading Jet $p_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'jeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='jeta', label='AK4 Leading Jet Eta'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'jphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='jphi', label='AK4 Leading Jet Phi'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'fjpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(38,250,2000, name='fjpt', label='AK15 Leading Jet $p_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'fjeta': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='fjeta', label='AK15 Leading Jet Eta'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'fjphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='fjphi', label='AK15 Leading Jet Phi'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'mT': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(60,0,1200, name='mT', label='$m_{T}$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'met': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(30,0,2000, name='met', label='MET $p_{T}$ [GeV]'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'metphi': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(35,-3.5,3.5, name='metphi', label='MET $\phi$'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				#hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'genWpt': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.Regular(50,0,1000, name='genWpt', label='gen W $p_{T}$'),
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
			'njets': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='njets', label='AK4 Number of Jets'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'njbtagL': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='njbtagL', label='AK4 Number of BTag Jets'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
			'nfjets': hist.Hist(
				hist.axis.StrCategory([], name='region', growth=True),
				hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='nfjets', label='AK15 Number of Jets'),
				hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
				storage=hist.storage.Weight(),
			),
	}

	def process(self, events):
		isData = not hasattr(events, "genWeight")
		if isData:
			# Nominal JEC are already applied in data
			return self.process_shift(events, None)
		
		# MC: nominal pass (all histograms; 'template' for "nominal" + every
		# weight variation) plus one pass per JES/JER shift ('template' only,
		# systematic=<shift>, nominal weight). The shifts change jet kinematics,
		# so the selection is rerun each time.
		return processor.accumulate(
			self.process_shift(events, shift) for shift in (None,) + tuple(self._jet_shifts)
		)
	


	def process_shift(self, events, shift_name):

		dataset = events.metadata['dataset']
		caller = dataset.split("____")[0]

		selected_regions = []
		for region, samples in self._samples.items():
			for sample in samples:
				if sample not in dataset: continue
				selected_regions.append(region)
		
		isData = not hasattr(events, "genWeight")
		# nominal JEC / nominal JER unless this pass is a JES or JER shift (MC only)
		jes_var, jer_var = self._jet_shifts.get(shift_name, ('nominal', 'nom'))
		doJER = (not isData) and (not self._skipJER)
		# Seed for the stochastic JER smearing: fixed per chunk, so every pass
		# (nominal, JES, JER) draws the same random numbers for the same jets.
		jer_seed = [int(events.run[0]), int(events.luminosityBlock[0]), int(events.event[0])] if len(events) > 0 else [0]
		def jer_factor(jets, pt, var, fat):
			"""JER factor for `jets` at JEC-corrected `pt`; 1 for data or _skipJER."""
			if not doJER:
				return ak.ones_like(pt)
			gen = jets.matched_gen
			genpt = ak.fill_none(gen.pt, 0.)
			gendr = ak.fill_none(jets.delta_r(gen), 99.)
			return get_jer_factor(self._year, pt, jets.eta, rho_density, genpt, gendr, var, fat, jer_seed + [15 if fat else 4])
		
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
		get_ele_trig_eff	  = self._corrections['get_ele_trig_eff']
		get_pho_trig_weight	  = self._corrections['get_pho_trig_weight']
		get_nlo_ewk_weight	   = self._corrections['get_nlo_ewk_weight']	
		#get_nnlo_nlo_weight	  = self._corrections['get_nnlo_nlo_weight'][self._year]
		get_btag_weight	  = self._corrections['get_btag_weight']
		get_ttag_weight	  = self._corrections['get_ttag_weight']
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
		## Per-electron HLT efficiency in data / MC (same Loose/Tight choice as hlt_sf),
		## for the 2-electron trigger SF. 0 below 25 GeV, where Ele30 cannot fire.
		eleEffD_loose, _, _, eleEffM_loose, _, _ = get_ele_trig_eff(self._year, abs(e.eta), e.pt, "Loose")
		eleEffD_tight, _, _, eleEffM_tight, _, _ = get_ele_trig_eff(self._year, abs(e.eta), e.pt, "Tight")
		e['hlt_eff_data'] = ak.where(
			e.istight,
			eleEffD_tight,
			ak.where(e.isloose, eleEffD_loose, ak.zeros_like(e.pt))
		)
		e['hlt_eff_mc'] = ak.where(
			e.istight,
			eleEffM_tight,
			ak.where(e.isloose, eleEffM_loose, ak.zeros_like(e.pt))
		)
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

		#resolvedT = events.TopResolved
		#ResolvedTvsQCD = resolvedT.TTScore / (resolvedT.TTScore + resolvedT.QCDScore)
		#resolvedT['ResolvedTvsQCD'] = resolvedT.TTScore / (resolvedT.TTScore + resolvedT.QCDScore)
		#nResolvedT = ak.num(resolvedT, axis=1)
		#leading_resolvedT = ak.firsts(resolvedT)

		fj = events.AK15PuppiJet
		rho_density = events.Rho.fixedGridRhoFastjetAll
		#fjec_corr = get_fjec_correction(self._year, fj.pt, fj.eta, fj.phi, rho_density, fj.area, run, isData)
		#fj['pt'] = fj.pt * fjec_corr
		#fj['mass'] = fj.mass * fjec_corr
		fj_pt_raw   = fj.pt   * (1 - fj.rawFactor)
		fj_mass_raw = fj.mass * (1 - fj.rawFactor)
		fjec_nom = get_fjec_correction(self._year, fj_pt_raw, fj.eta, fj.phi, rho_density, fj.area, run, isData)
		fjec_corr = fjec_nom
		if jes_var != 'nominal':
			fjec_corr = get_fjec_correction(self._year, fj_pt_raw, fj.eta, fj.phi, rho_density, fj.area, run, isData, var=jes_var)
		# JER smearing factor, always derived from the nominal-JEC pt (so a JES
		# pass uses the same smearing as the nominal pass)
		fjer = jer_factor(fj, fj_pt_raw * fjec_nom, jer_var, True)
		fj['pt']   = fj_pt_raw   * fjec_corr * fjer
		fj['mass'] = fj_mass_raw * fjec_corr * fjer
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
		#jec_corr = get_jec_correction(self._year, j.pt, j.eta, j.phi, rho_density, j.area, run, isData)
		#j['pt'] = j.pt * jec_corr
		#j['mass'] = j.mass * jec_corr
		j_pt_raw   = j.pt   * (1 - j.rawFactor)
		j_mass_raw = j.mass * (1 - j.rawFactor)
		jec_nom = get_jec_correction(self._year, j_pt_raw, j.eta, j.phi, rho_density, j.area, run, isData)
		jec_corr = jec_nom
		if jes_var != 'nominal':
			jec_corr = get_jec_correction(self._year, j_pt_raw, j.eta, j.phi, rho_density, j.area, run, isData, var=jes_var)
		jer = jer_factor(j, j_pt_raw * jec_nom, jer_var, False)
		if shift_name in self._jet_shifts:
			# Propagate the AK4 JES/JER shift to MET (simplified type-1: every AK4
			# jet with nominal pt > 15 GeV): MET -= sum(shifted - nominal), where
			# nominal = nominal JEC x nominal JER. The nominal MET itself is not
			# re-derived (MET_TODO.txt). met is only used after this point.
			j_pt_nom = j_pt_raw * jec_nom * (jer if jer_var == 'nom' else jer_factor(j, j_pt_raw * jec_nom, 'nom', False))
			dpt = ak.where(j_pt_nom > 15., j_pt_raw * jec_corr * jer - j_pt_nom, 0.)
			met_x = met.pt * np.cos(met.phi) - ak.sum(dpt * np.cos(j.phi), axis=1)
			met_y = met.pt * np.sin(met.phi) - ak.sum(dpt * np.sin(j.phi), axis=1)
			met['pt'], met['phi'] = np.hypot(met_x, met_y), np.arctan2(met_y, met_x)
		j['pt']   = j_pt_raw   * jec_corr * jer
		j['mass'] = j_mass_raw * jec_corr * jer
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
			
			gen = events.GenPart

			###
			# AK15 truth class for top-tag weight (top / b / qcd), from topstudy_run.py
			###
			#j['isbflav'] = j.hadronFlavour == 5
			#j_bflav = j[j.isbflav]
			#binfj = leading_fj.vec.delta_r(j_bflav)
			#n_binfj = ak.sum(binfj < 1.5, axis=1)
			n_binfj = ak.fill_none(leading_fj.nBHadrons, 0)

			#yourmotheris = gen.distinctParent.distinctParent.pdgId
			gen['isq'] = (abs(gen.pdgId) >= 1) & (abs(gen.pdgId) <= 5) & gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			p1 = ak.fill_none(gen.distinctParent.pdgId, 0)
			#p2 = ak.fill_none(ak.fill_none(gen.distinctParent.distinctParent.pdgId, 0), 0)
			p2 = ak.fill_none(gen.distinctParent.distinctParent.pdgId, 0)
			gen['isqfromt'] = gen.isq & (abs(p1) == 24) & (abs(p2) == 6)
			gen['isbfromt'] = gen.isq & (abs(gen.pdgId) == 5) & (abs(p1) == 6)
			q_from_t = gen[gen.isqfromt]
			b_from_t = gen[gen.isbfromt]

			# keep only b's whose top also decayed hadronically (2 W-quarks share its index)
			b_top = b_from_t.distinctParentIdxG
			q_top = q_from_t.distinctParent.distinctParentIdxG
			b_from_t = b_from_t[ak.sum(b_top[:, :, None] == q_top[:, None, :], axis=2) >= 2]
			q_from_t = ak.concatenate([q_from_t, b_from_t], axis=1)

			qinfj = leading_fj.vec.delta_r(q_from_t)
			n_qinfj = ak.sum(qinfj <= 1.5, axis=1)
			#fj['qisinfj'] = ak.all(fj.vec.metric_table(q_from_t)<1.5, axis=2)
			#q_from_t_fj = fj[fj.qisinfj]
			#nq_from_t_fj = ak.num(q_from_t_fj, axis=1)

			# None if no AK15 jet
			fj_class = ak.where(n_qinfj >= 2, "top", ak.where(n_binfj > 0, "b", "qcd"))

			gen['isb'] =(abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
			genTops = gen[gen.isTop]
			nlo = np.ones(len(events), dtype='float')
			if('TTto' in dataset): 
				nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt) * get_ttbar_weight(genTops[:,1].pt))
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
			GenIsoPho = events.GenIsolatedPhoton
			nGenIsoPho = ak.num(GenIsoPho, axis=1)
			#
			nlo_ewk = np.ones(len(events), dtype='float')
			if('WtoLNu' in dataset): 
				nlo_ewk = get_nlo_ewk_weight['w'](ak.max(genWs.pt, axis=1))
			elif('DY' in dataset): 
				nlo_ewk = get_nlo_ewk_weight['dy'](ak.max(genDYs.pt, axis=1))
			elif('Zto2Nu' in dataset): 
				nlo_ewk = get_nlo_ewk_weight['z'](ak.max(genZs.pt, axis=1))
			elif('GJ' in dataset):
				nlo_ewk = get_nlo_ewk_weight['a'](ak.max(genAs.pt, axis=1))

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
			# AK15 top-tagging weights
			###
			ttagSF, ttagSFUp, ttagSFDown = get_ttag_weight(self._year).ttag_weight(
				leading_fj.pt,
				fj_class,
				leading_fj.TvsQCD > self._TvsQCDwp[self._year]
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
			weights.add('genw', events.genWeight)
			# Weights.add(name, w, wUp, wDown): wUp/wDown are absolute (coffea
			# divides by w where w != 0); variations appear as <name>Up/Down.
			# Same sources and names as processors/topstudy_unc.py.
			nev = len(events)
			one = np.ones(nev)
			pu_nom = finite_or(pu, 1.)
			weights.add('pileup', pu_nom, finite_or(pu_up, pu_nom), finite_or(pu_down, pu_nom))
			weights.add('nlo_ewk',nlo_ewk)
			# top-pT reweighting: up = applied twice, down = not applied
			# (nlo = 1 for non-ttbar samples, so no effect there)
			nlo_nom = finite_or(ak.to_numpy(nlo) if not isinstance(nlo, np.ndarray) else nlo, 1.)
			weights.add('topptreweighting', nlo_nom, nlo_nom**2, one)
			# parton shower, muR/muF scale envelope, PDF and alpha_S (relative
			# weights; 1 where the sample has no such branch)
			(isr_up, isr_down), (fsr_up, fsr_down) = get_ps_weights(events)
			weights.add('ISR', one, isr_up, isr_down)
			weights.add('FSR', one, fsr_up, fsr_down)
			scale_up, scale_down = get_scale_weights(events)
			weights.add('QCDscale', one, scale_up, scale_down)
			(pdf_up, pdf_down), (as_up, as_down) = get_pdf_weights(events)
			weights.add('PDF', one, pdf_up, pdf_down)
			weights.add('alphaS', one, as_up, as_down)
			weights.add('btagSF',btagSF)
			# b-tag variations as ratios to the nominal btagSF. btag_weight()
			# returns btagSF = 0 when an untagged jet has SF*eff >= 1; a plain
			# Up/btagSF is then 0/0 = NaN and poisons whole bins of every
			# variation. Those events have zero nominal weight, so their
			# variation ratio is set to 1.
			btag_nom = np.asarray(ak.to_numpy(btagSF), dtype=np.float64)
			btag_ok = np.isfinite(btag_nom) & (btag_nom > 0)
			btag_den = np.where(btag_ok, btag_nom, 1.)
			def btag_ratio(var):
				r = np.asarray(ak.to_numpy(var), dtype=np.float64) / btag_den
				return np.where(btag_ok & np.isfinite(r), r, 1.)
			weights.add('btagSFbc_correlated', one, btag_ratio(btagSFbc_correlatedUp), btag_ratio(btagSFbc_correlatedDown))
			weights.add('btagSFbc_uncorrelated', one, btag_ratio(btagSFbc_uncorrelatedUp), btag_ratio(btagSFbc_uncorrelatedDown))
			weights.add('btagSFlight_correlated', one, btag_ratio(btagSFlight_correlatedUp), btag_ratio(btagSFlight_correlatedDown))
			weights.add('btagSFlight_uncorrelated', one, btag_ratio(btagSFlight_uncorrelatedUp), btag_ratio(btagSFlight_uncorrelatedDown))
			weights.add('ttagSF', ttagSF, ttagSFUp, ttagSFDown)

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
			# photon SFs vary independently (masks are disjoint). The nominal
			# product of these weights equals the former 'ids'*'reco'*'iso'*'trig'.
			# As before, the 2-lepton regions use the product over the TIGHT
			# leptons (the nominal is kept unchanged here; topstudy_unc.py uses
			# the loose leptons there).
			###
			def obj_prod(objs, field):
				return finite_or(ak.to_numpy(ak.fill_none(ak.prod(objs[field], axis=1), 1.0)), 1.)
			def add_obj_sf(name, pieces):
				"""pieces: list of (mask, objects, field). Builds nominal/up/down
				from field, field+'_up', field+'_down' on the masked events."""
				w = {sfx: np.ones(nev) for sfx in ('', '_up', '_down')}
				for mask, objs, field in pieces:
					for sfx in w:
						w[sfx][mask] = obj_prod(objs, field+sfx)[mask]
				weights.add(name, w[''], w['_up'], w['_down'])

			# ID: electron, muon, photon
			add_obj_sf('ids_e',   [(mask_isoneE, e_tight, 'id_sf'), (mask_istwoE, e_tight, 'id_sf')])
			add_obj_sf('ids_mu',  [(mask_isoneM, mu_tight, 'id_sf'), (mask_istwoM, mu_tight, 'id_sf')])
			add_obj_sf('ids_pho', [(mask_isoneG, pho_tight, 'id_sf')])
			# RECO: electron
			add_obj_sf('reco_e',  [(mask_isoneE, e_tight, 'reco_sf'), (mask_istwoE, e_tight, 'reco_sf')])
			# ISO: muon
			add_obj_sf('iso_mu',  [(mask_isoneM, mu_tight, 'iso_sf'), (mask_istwoM, mu_tight, 'iso_sf')])

			###
			# Trigger efficiency weights: electron (1e, 2e), MET (0l, 1mu, 2mu), photon (1g)
			###
			met_sf_0,  met_sf_0_up,  met_sf_0_down  = get_met_trig_weight(self._year, u['sr'].r)
			met_sf_m,  met_sf_m_up,  met_sf_m_down  = get_met_trig_weight(self._year, u['wmcr'].r)
			met_sf_mm, met_sf_mm_up, met_sf_mm_down = get_met_trig_weight(self._year, u['zmcr'].r)
			pho_sf, pho_sf_up, pho_sf_down = get_pho_trig_weight(self._year, ak.firsts(pho_tight).pt)
			def ev_arr(x):
				return finite_or(ak.to_numpy(ak.fill_none(x, 1.0)), 1.)

			## 2e: the single-e trigger fires if EITHER electron fires it:
			## SF = [1 - prod(1 - eff_data)] / [1 - prod(1 - eff_mc)]
			trg_p_data = 1 - ak.prod(1 - e_loose.hlt_eff_data, axis=1)
			trg_p_mc   = 1 - ak.prod(1 - e_loose.hlt_eff_mc,   axis=1)
			trg_sf_istwoE = ev_arr(ak.where(trg_p_mc > 0, trg_p_data / trg_p_mc, 1.0))
			## 2e up/down: nominal x R, where R is the change of the same OR
			## formula when each electron's data efficiency eff_mc*SF follows its
			## per-electron SF up/down (the uncertainty the 1e region uses).
			def trg_p_sf(field):
				return 1 - ak.prod(1 - np.minimum(e_loose.hlt_eff_mc * e_loose[field], 1.), axis=1)
			trg_p_nom = trg_p_sf('hlt_sf')
			def trg_2e_ratio(field):
				return ev_arr(ak.where(trg_p_nom > 0, trg_p_sf(field) / trg_p_nom, 1.0))

			trg_e = {sfx: np.ones(nev) for sfx in ('', '_up', '_down')}
			for sfx in trg_e:
				trg_e[sfx][mask_isoneE] = obj_prod(e_tight, 'hlt_sf'+sfx)[mask_isoneE]
				twoE = trg_sf_istwoE if sfx == '' else trg_sf_istwoE * trg_2e_ratio('hlt_sf'+sfx)
				trg_e[sfx][mask_istwoE] = twoE[mask_istwoE]
			weights.add('trig_e', trg_e[''], trg_e['_up'], trg_e['_down'])

			trg_met = {sfx: np.ones(nev) for sfx in ('', '_up', '_down')}
			for mask, trio in ((mask_iszeroL, (met_sf_0,  met_sf_0_up,  met_sf_0_down)),
			                   (mask_isoneM,  (met_sf_m,  met_sf_m_up,  met_sf_m_down)),
			                   (mask_istwoM,  (met_sf_mm, met_sf_mm_up, met_sf_mm_down))):
				for sfx, x in zip(('', '_up', '_down'), trio):
					trg_met[sfx][mask] = ev_arr(x)[mask]
			weights.add('trig_met', trg_met[''], trg_met['_up'], trg_met['_down'])

			trg_pho = {sfx: np.ones(nev) for sfx in ('', '_up', '_down')}
			for sfx, x in zip(('', '_up', '_down'), (pho_sf, pho_sf_up, pho_sf_down)):
				trg_pho[sfx][mask_isoneG] = ev_arr(x)[mask_isoneG]
			weights.add('trig_pho', trg_pho[''], trg_pho['_up'], trg_pho['_down'])



		def normalize(val, cut,fillvalue=np.nan):
			if cut is None:
				ar = ak.to_numpy(ak.fill_none(val, fillvalue))
				return ar
			else:
				ar = ak.to_numpy(ak.fill_none(val[cut], fillvalue))
				return ar
		#def normalize(val, cut):
		#	if cut is None:
		#		ar = ak.to_numpy(ak.fill_none(val, np.nan))
		#		return ar
		#	else:
		#		ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
		#		return ar
				
		def fill(region, systematic):
			cut = selection.all(*regions[region])
			sname = 'nominal' if systematic is None else systematic
			if systematic in weights.variations:
				weight = weights.weight(modifier=systematic)[cut]
			else:
				weight = weights.weight()[cut]
			#print('weights: ', weights.weightStatistics.keys())
			# filled for every systematic (the histograms below are nominal only)
			output['template'].fill(
				  region=region,
				  systematic=sname,
				  recoil=normalize(u[region].r, cut),
				  TvsQCD=normalize(leading_fj.TvsQCD, cut),
				  weight=weight
			)
			if systematic is None:
				variables = {
					'met':					met.pt,
					'metphi':				met.phi,
					#'ut' :					u[region].r,
					'uphi' :				u[region].phi,
					'mindphirecoil':		ak.min(abs(u['sr'].delta_phi(j_clean.T)), axis=1,mask_identity=False),
					'minDphirecoil':		abs(u[region].delta_phi(leading_fj.T)),
					'mupt':					leading_mu.pt,
					'muphi':				leading_mu.phi,
					'mueta':				leading_mu.eta,
					'elept':				leading_e.pt,
					'elephi':				leading_e.phi,
					'eleeta':				leading_e.eta,
					'phopt':				leading_pho.pt,
					'phophi':				leading_pho.phi,
					'phoeta':				leading_pho.eta,
					'fjpt':					leading_fj.pt,
					'fjphi':				leading_fj.phi,
					'fjeta':				leading_fj.eta,
					'jpt':					leading_j.pt,
					'jphi':					leading_j.phi,
					'jeta':					leading_j.eta,
					'njbtagL':				j_nbtagvL,
					'njets':				j_nclean,
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
						variables['dimuphi']	= dimu.phi
						variables['dimueta']	= dimu.eta

				for variable in output:
					if variable not in variables:
						continue
					normalized_variable = {variable: normalize(variables[variable],cut)}
					output[variable].fill(
						region=region,
						TvsQCD=normalize(leading_fj.TvsQCD, cut),
						**normalized_variable,
						weight=weight,
					)
				if not isData:
					output['genWpt'].fill(
						region=region,
						genWpt=normalize(ak.firsts(genW.pt), cut),
						weight=weight
					)
					output['fjpt_cat'].fill(
						  region=region,
						  fjpt=normalize(leading_fj.pt, cut),
						  TvsQCD=normalize(leading_fj.TvsQCD, cut),
						  n_qinfj=normalize(n_qinfj, cut,fillvalue=-1),
						  n_binfj=normalize(n_binfj, cut,fillvalue=-1),
						  weight=weight
					)
				output['ut'].fill(
					  region=region,
					  ut=normalize(u[region].r, cut),
					  TvsQCD=normalize(leading_fj.TvsQCD, cut),
					  #ResolvedTvsQCD=normalize(leading_resolvedT.ResolvedTvsQCD, cut),
					  weight=weight
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
			## Store cutflow (nominal pass only; JES/JER passes fill 'template' only)
			#########
			if shift_name is not None:
				continue
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
	# `python3 processors/dev_run3.py` puts processors/ (not analysis/) on
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
