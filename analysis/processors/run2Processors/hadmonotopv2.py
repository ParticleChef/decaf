#!/usr/bin/env python
import logging
import numpy as np
import awkward as ak
import json
import copy
from collections import defaultdict
from coffea import processor
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
        #Values from https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2    
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVDatasetsUL2016
        '2016postVFP': 16.8,
        '2016preVFP': 19.5,
        '2017': 41.48,
        '2018': 59.83
    }

    lumiMasks = {
        '2016postVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2016preVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2017': LumiMask("data/jsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
        '2018"': LumiMask("data/jsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
    }
    
    met_filters = {
        # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
        '2016postVFP': [
                'goodVertices',
                'globalSuperTightHalo2016Filter',
                'HBHENoiseFilter',
                'HBHENoiseIsoFilter',
                'EcalDeadCellTriggerPrimitiveFilter',
                'BadPFMuonFilter',
                'BadPFMuonDzFilter',
                'eeBadScFilter'
                ],

        '2016preVFP': [
                'goodVertices',
                'globalSuperTightHalo2016Filter',
                'HBHENoiseFilter',
                'HBHENoiseIsoFilter',
                'EcalDeadCellTriggerPrimitiveFilter',
                'BadPFMuonFilter',
                'BadPFMuonDzFilter',
                'eeBadScFilter'
                ],
        
        '2017': [
                'goodVertices', 
                'globalSuperTightHalo2016Filter', 
                'HBHENoiseFilter', 
                'HBHENoiseIsoFilter', 
                'EcalDeadCellTriggerPrimitiveFilter', 
                'BadPFMuonFilter', 
                'BadPFMuonDzFilter', 
                'eeBadScFilter', 
                'ecalBadCalibFilter'
                ],

        '2018': [
                'goodVertices', 
                'globalSuperTightHalo2016Filter', 
                'HBHENoiseFilter', 
                'HBHENoiseIsoFilter', 
                'EcalDeadCellTriggerPrimitiveFilter', 
                'BadPFMuonFilter', 
                'BadPFMuonDzFilter', 
                'eeBadScFilter', 
                'ecalBadCalibFilter'
                ]
    }
            
    def __init__(self, year, xsec, corrections, ids, common):

        self._year = year
        self._lumi = 1000.*float(AnalysisProcessor.lumis[year])
        self._xsec = xsec
        self._systematics = True
        self._skipJER = False

        self._samples = {
            'sr':('Z1Jets','Z2Jets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET','TPhiTo2Chi'),
            'wmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'tmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'wecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'tecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'qcdcr':('Z1Jets','Z2Jets','WJets','TT','ST','WW','WZ','ZZ','QCD','MET'),
        }
        
        self._TvsQCDwp = {
            '2016preVFP': 0.53,
            '2016postVFP': 0.53,
            '2017': 0.61,
            '2018': 0.65
        }

        self._met_triggers = {
            '2016postVFP': [
                'PFMETNoMu90_PFMHTNoMu90_IDTight',
                'PFMETNoMu100_PFMHTNoMu100_IDTight',
                'PFMETNoMu110_PFMHTNoMu110_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2016preVFP': [
                'PFMETNoMu90_PFMHTNoMu90_IDTight',
                'PFMETNoMu100_PFMHTNoMu100_IDTight',
                'PFMETNoMu110_PFMHTNoMu110_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2017': [
                'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2018': [
                'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ]
        }

        self._singleelectron_triggers = { #2017 and 2018 from monojet, applying dedicated trigger weights
            '2016postVFP': [
                'Ele27_WPTight_Gsf',
                'Ele105_CaloIdVT_GsfTrkIdT'
            ],
            '2016preVFP': [
                'Ele27_WPTight_Gsf',
                'Ele105_CaloIdVT_GsfTrkIdT'
            ],
            '2017': [
                'Ele35_WPTight_Gsf',
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Photon200'
            ],
            '2018': [
                'Ele32_WPTight_Gsf',
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Photon200'
            ]
        }

        self._corrections = corrections
        self._ids = ids
        self._common = common

        ptbins=[30.0, 
                60.0, 
                90.0, 
                120.0, 
                150.0, 
                180.0, 
                210.0, 
                250.0, 
                280.0, 
                310.0, 
                340.0, 
                370.0, 
                400.0, 
                430.0, 
                470.0, 
                510.0, 
                550.0, 
                590.0, 
                640.0, 
                690.0, 
                740.0, 
                790.0, 
                840.0, 
                900.0, 
                960.0, 
                1020.0, 
                1090.0, 
                1160.0, 
                1250.0]

        self.make_output = lambda: {
            'sumw': 0.,
            'template': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.StrCategory([], name='systematic', growth=True),
                hist.axis.Variable([250,310,370,470,590,840,1020,1250,3000], name='recoil', label=r'$U$ [GeV]'),
                hist.axis.Variable([40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300], name='fjmass', label=r'AK15 Jet $m_{sd}$'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'TvsQCD': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(15,0,1, name='TvsQCD', label='TvsQCD'),
                storage=hist.storage.Weight(),
            ),
            'mindphirecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='mindphirecoil', label='Min dPhi(Recoil,AK4s)'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'minDphirecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='minDphirecoil', label='Min dPhi(Recoil,AK15s)'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'CaloMinusPfOverRecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,0,1, name='CaloMinusPfOverRecoil', label='Calo - Pf / Recoil'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'met': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,600, name='met', label='MET'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'metphi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='metphi', label='MET phi'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'mindphimet': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='mindphimet', label='Min dPhi(MET,AK4s)'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'minDphimet': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='minDphimet', label='Min dPhi(MET,AK15s)'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'j1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Variable(ptbins, name='j1pt', label='AK4 Leading Jet Pt'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'j1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1eta', label='AK4 Leading Jet Eta'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'j1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1phi', label='AK4 Leading Jet Phi'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'fj1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Variable(ptbins, name='fj1pt', label='AK15 Leading SoftDrop Jet Pt'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'fj1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1eta', label='AK15 Leading SoftDrop Jet Eta'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'fj1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1phi', label='AK15 Leading SoftDrop Jet Phi'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'njets': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='njets', label='AK4 Number of Jets'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'ndflvL': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='ndflvL', label='AK4 Number of deepFlavor Loose Jets'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'nfjclean': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4], name='nfjclean', label='AK15 Number of Cleaned Jets'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'mT': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(20,0,600, name='mT', label='Transverse Mass'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'l1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Variable(ptbins, name='l1pt', label='Leading Lepton Pt'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'l1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(48,-2.4,2.4, name='l1eta', label='Leading Lepton Eta'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'l1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(64,-3.2,3.2, name='l1phi', label='Leading Lepton Phi'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
    }

    def process(self, events):
        isData = not hasattr(events, "genWeight")
        if isData:
            # Nominal JEC are already applied in data
            return self.process_shift(events, None)

        jet_factory              = self._corrections['jet_factory']
        #fatjet_factory           = self._corrections['fatjet_factory']
        subjet_factory           = self._corrections['subjet_factory']
        met_factory              = self._corrections['met_factory']

        import cachetools
        jec_cache = cachetools.Cache(np.inf)
    
        nojer = "NOJER" if self._skipJER else ""
        thekey = f"{self._year}mc{nojer}"

        def add_jec_variables(jets, event_rho):
            jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
            jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
            jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
            jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
            return jets
        
        jets = jet_factory[thekey].build(add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll), jec_cache)
        #fatjets = fatjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiJet, events.fixedGridRhoFastjetAll), jec_cache)
        subjets = subjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiSubJet, events.fixedGridRhoFastjetAll), jec_cache)
        met = met_factory.build(events.MET, jets, {})

        shifts = [({"Jet": jets, "AK15PFPuppiSubJet": subjets, "MET": met}, None)]
        if self._systematics:
            shifts.extend([
                ({"Jet": jets.JES_jes.up, "AK15PFPuppiSubJet": subjets.JES_jes.up, "MET": met.JES_jes.up}, "JESUp"),
                ({"Jet": jets.JES_jes.down, "AK15PFPuppiSubJet": subjets.JES_jes.down, "MET": met.JES_jes.down}, "JESDown"),
                ({"Jet": jets, "AK15PFPuppiSubJet": subjets, "MET": met.MET_UnclusteredEnergy.up}, "UESUp"),
                ({"Jet": jets, "AK15PFPuppiSubJet": subjets, "MET": met.MET_UnclusteredEnergy.down}, "UESDown"),
            ])
            if not self._skipJER:
                shifts.extend([
                    ({"Jet": jets.JER.up, "AK15PFPuppiSubJet": subjets.JER.up, "MET": met.JER.up}, "JERUp"),
                    ({"Jet": jets.JER.down, "AK15PFPuppiSubJet": subjets.JER.down, "MET": met.JER.down}, "JERDown"),
                ])
        return processor.accumulate(self.process_shift(update(events, collections), name) for collections, name in shifts)

    def process_shift(self, events, shift_name):

        dataset = events.metadata['dataset']

        selected_regions = []
        for region, samples in self._samples.items():
            for sample in samples:
                if sample not in dataset: continue
                selected_regions.append(region)

        isData = not hasattr(events, "genWeight")
        selection = PackedSelection(dtype="uint64")
        weights = Weights(len(events), storeIndividual=True)
        output = self.make_output()
        if shift_name is None and not isData:
            output['sumw'] = ak.sum(events.genWeight)

        ###
        #Getting corrections, ids from .coffea files
        ###

        get_met_trig_weight      = self._corrections['get_met_trig_weight']
        get_ele_loose_id_sf      = self._corrections['get_ele_loose_id_sf']
        get_ele_tight_id_sf      = self._corrections['get_ele_tight_id_sf']
        get_ele_trig_weight      = self._corrections['get_ele_trig_weight']
        get_ele_reco_sf_below20  = self._corrections['get_ele_reco_sf_below20']
        get_ele_reco_sf_above20  = self._corrections['get_ele_reco_sf_above20']
        get_mu_loose_id_sf       = self._corrections['get_mu_loose_id_sf']
        get_mu_tight_id_sf       = self._corrections['get_mu_tight_id_sf']
        get_mu_loose_iso_sf      = self._corrections['get_mu_loose_iso_sf']
        get_mu_tight_iso_sf      = self._corrections['get_mu_tight_iso_sf']
        get_mu_rochester_sf      = self._corrections['get_mu_rochester_sf'][self._year]
        get_met_xy_correction    = self._corrections['get_met_xy_correction']
        get_pu_weight            = self._corrections['get_pu_weight']    
        get_nlo_ewk_weight       = self._corrections['get_nlo_ewk_weight']    
        get_nnlo_nlo_weight      = self._corrections['get_nnlo_nlo_weight'][self._year]
        get_msd_corr             = self._corrections['get_msd_corr']
        get_btag_weight      = self._corrections['get_btag_weight']
        
        isLooseElectron = self._ids['isLooseElectron'] 
        isTightElectron = self._ids['isTightElectron'] 
        isLooseMuon     = self._ids['isLooseMuon']     
        isTightMuon     = self._ids['isTightMuon']     
        isLooseTau      = self._ids['isLooseTau']      
        isLoosePhoton   = self._ids['isLoosePhoton']   
        isTightPhoton   = self._ids['isTightPhoton']   
        isGoodAK4       = self._ids['isGoodAK4']       
        isGoodAK15    = self._ids['isGoodAK15']    
        isHEMJet        = self._ids['isHEMJet']        
        
        deepflavWPs = self._common['btagWPs']['deepflav'][self._year]
        deepcsvWPs = self._common['btagWPs']['deepcsv'][self._year]

        ###
        #Initialize global quantities (MET ecc.)
        ###

        npv = events.PV.npvsGood
        run = events.run
        calomet = events.CaloMET
        met = events.MET
        met['pt'] , met['phi'] = get_met_xy_correction(self._year, npv, run, met.pt, met.phi, isData)

        ###
        #Initialize physics objects
        ###

        mu = events.Muon
        if isData:
            k = get_mu_rochester_sf.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi)
        else:
            kspread = get_mu_rochester_sf.kSpreadMC(
                mu.charge, 
                mu.pt, 
                mu.eta, 
                mu.phi,
                mu.matched_gen.pt
            )
            mc_rand = ak.unflatten(np.random.rand(ak.count(ak.flatten(mu.pt))), ak.num(mu.pt))
            ksmear = get_mu_rochester_sf.kSmearMC(
                mu.charge, 
                mu.pt, 
                mu.eta, 
                mu.phi,
                mu.nTrackerLayers, 
                mc_rand
            )
            hasgen = ~np.isnan(ak.fill_none(events.Muon.matched_gen.pt, np.nan))
            k = ak.where(hasgen, kspread, ksmear)
        mu['pt'] = ak.where((mu.pt<200),k*mu.pt, mu.pt)
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
            mu.iso_sf
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
        leading_mu = ak.firsts(mu_tight)

        e = events.Electron
        e['isclean'] = ak.all(e.metric_table(mu_loose) > 0.3, axis=2)
        e['reco_sf'] = ak.where(
            (e.pt<20),
            get_ele_reco_sf_below20(self._year, e.eta+e.deltaEtaSC, e.pt), 
            get_ele_reco_sf_above20(self._year, e.eta+e.deltaEtaSC, e.pt)
        )
        e['isloose'] = isLooseElectron(e,self._year)
        e['id_sf'] = ak.where(
            e.isloose,
            get_ele_loose_id_sf(self._year, e.eta+e.deltaEtaSC, e.pt),
            ak.ones_like(e.pt)
        )
        e['istight'] = isTightElectron(e,self._year)
        e['id_sf'] = ak.where(
            e.istight,
            get_ele_tight_id_sf(self._year, e.eta+e.deltaEtaSC, e.pt),
            e.id_sf
        )
        e['T'] = ak.zip(
            {
                "r": e.pt,
                "phi": e.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        e_clean = e[e.isclean]
        e_loose = e_clean[e_clean.isloose]
        e_tight = e_clean[e_clean.istight]
        e_ntot = ak.num(e, axis=1)
        e_nloose = ak.num(e_loose, axis=1)
        e_ntight = ak.num(e_tight, axis=1)
        leading_e = ak.firsts(e_tight)

        tau = events.Tau
        tau['isclean']=(
            ak.all(tau.metric_table(mu_loose) > 0.4, axis=2) 
            & ak.all(tau.metric_table(e_loose) > 0.4, axis=2)
        )
        
        tau['isloose']=isLooseTau(tau, self._year)
        tau_clean=tau[tau.isclean]
        tau_loose=tau_clean[tau_clean.isloose]
        tau_ntot=ak.num(tau, axis=1)
        tau_nloose=ak.num(tau_loose, axis=1)

        pho = events.Photon
        pho['isclean']=(
            ak.all(pho.metric_table(mu_loose) > 0.5, axis=2)
            & ak.all(pho.metric_table(e_loose) > 0.5, axis=2)
            & ak.all(pho.metric_table(tau_loose) > 0.5, axis=2)
        )
        pho['isloose']=isLoosePhoton(pho,self._year)
        pho_clean=pho[pho.isclean]
        pho_loose=pho_clean[pho_clean.isloose]
        pho_ntot=ak.num(pho, axis=1)
        pho_nloose=ak.num(pho_loose, axis=1)

        fj = events.AK15PFPuppiJet
        fj['pt'] = fj.subjets.sum().pt
        fj['msd_corr'] = get_msd_corr(fj)
        fj['isclean'] = (
            ak.all(fj.metric_table(mu_loose) > 1.5, axis=2)
            & ak.all(fj.metric_table(e_loose) > 1.5, axis=2)
            & ak.all(fj.metric_table(tau_loose) > 1.5, axis=2)
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
        probQCD=fj.particleNetAK15_QCDbb+fj.particleNetAK15_QCDcc+fj.particleNetAK15_QCDb+fj.particleNetAK15_QCDc+fj.particleNetAK15_QCDothers
        probT=fj.particleNetAK15_Tbqq+fj.particleNetAK15_Tbcq
        fj['TvsQCD'] = probT/(probT+probQCD)
        fj_good = fj[fj.isgood]
        fj_clean = fj_good[fj_good.isclean]
        fj_ntot = ak.num(fj, axis=1)
        fj_ngood = ak.num(fj_good, axis=1)
        fj_nclean = ak.num(fj_clean, axis=1)
        leading_fj = ak.firsts(fj_clean)

        j = events.Jet
        j['isgood'] = isGoodAK4(j, self._year)
        j['isHEM'] = isHEMJet(j)
        j['isclean'] = (
            ak.all(j.metric_table(mu_loose) > 0.4, axis=2)
            & ak.all(j.metric_table(e_loose) > 0.4, axis=2)
            & ak.all(j.metric_table(tau_loose) > 0.4, axis=2)
            & ak.all(j.metric_table(pho_loose) > 0.4, axis=2)
        )
        j['isiso'] = ak.all(j.metric_table(leading_fj) > 1.5, axis=2)
        j['isdcsvL'] = (j.btagDeepB>deepcsvWPs['loose'])
        j['isdflvL'] = (j.btagDeepFlavB>deepflavWPs['loose'])
        j['T'] = ak.zip(
            {
                "r": j.pt,
                "phi": j.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        j_good = j[j.isgood]
        j_clean = j_good[j_good.isclean]
        j_iso = j_clean[j_clean.isiso]
        j_dcsvL = j_iso[j_iso.isdcsvL]
        j_dflvL = j_iso[j_iso.isdflvL]
        j_HEM = j[j.isHEM]
        j_ntot=ak.num(j, axis=1)
        j_ngood=ak.num(j_good, axis=1)
        j_nclean=ak.num(j_clean, axis=1)
        j_niso=ak.num(j_iso, axis=1)
        j_ndcsvL=ak.num(j_dcsvL, axis=1)
        j_ndflvL=ak.num(j_dflvL, axis=1)
        j_nHEM = ak.num(j_HEM, axis=1)
        leading_j = ak.firsts(j_clean)

        ###
        # Calculate recoil and transverse mass
        ###

        u = {
            'sr'    : met,
            'wecr'  : met+leading_e.T,
            'tecr'  : met+leading_e.T,
            'wmcr'  : met+leading_mu.T,
            'tmcr'  : met+leading_mu.T,
            'qcdcr' : met,
        }

        mT = {
            'wecr'  : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            'tecr'  : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            'wmcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
            'tmcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))) 
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
            nlo = np.ones(len(events), dtype='float')
            if('TTJets' in dataset): 
                nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt) * get_ttbar_weight(genTops[:,1].pt))
                
            gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            
            genWs = gen[gen.isW] 
            genZs = gen[gen.isZ]
            genDYs = gen[gen.isZ&(gen.mass>30)]
            
            nnlo_nlo = {}
            nlo_qcd = np.ones(len(events), dtype='float')
            nlo_ewk = np.ones(len(events), dtype='float')
            if('WJets' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['w'](genWs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['w'](genWs.pt.max())
                for systematic in get_nnlo_nlo_weight['w']:
                    nnlo_nlo[systematic] = ak.where(
                        ((ak.num(genWs, axis=1)>0)&(ak.firsts(genWs).pt>=100)),
                        get_nnlo_nlo_weight['w'][systematic](ak.firsts(genWs).pt),
                        np.ones(len(events), dtype='float')
                    )
            elif('DY' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['dy'](genDYs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['dy'](genDYs.pt.max())
                for systematic in get_nnlo_nlo_weight['dy']:
                    nnlo_nlo[systematic] = ak.where(
                        ((ak.num(genDYs, axis=1)>0)&(ak.firsts(genDYs).pt>=100)),
                        get_nnlo_nlo_weight['dy'][systematic](ak.firsts(genDYs).pt),
                        np.ones(len(events), dtype='float')
                    )
            elif('Z1Jets' in dataset or 'Z2Jets' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['z'](genZs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['z'](genZs.pt.max())
                for systematic in get_nnlo_nlo_weight['z']:
                    nnlo_nlo[systematic] = ak.where(
                        ((ak.num(genZs, axis=1)>0)&(ak.firsts(genZs).pt>=100)),
                        get_nnlo_nlo_weight['z'][systematic](ak.firsts(genZs).pt),
                        np.ones(len(events), dtype='float')
                    )
            ###
            # Calculate PU weight and systematic variations
            ###

            pu = get_pu_weight(self._year, events.Pileup.nTrueInt)

            ###
            # Trigger efficiency weight
            ###
           
            trig = {
                'sr':   get_met_trig_weight(self._year, met.pt),
                'wmcr': get_met_trig_weight(self._year, u['wmcr'].r),
                'tmcr': get_met_trig_weight(self._year, u['tmcr'].r),
                'wecr': get_ele_trig_weight(self._year, leading_e.eta+leading_e.deltaEtaSC, leading_e.pt),
                'tecr': get_ele_trig_weight(self._year, leading_e.eta+leading_e.deltaEtaSC, leading_e.pt),
                'qcdcr': get_met_trig_weight(self._year, met.pt),
            }

            ### 
            # Calculating electron and muon ID weights
            ###

            ids ={
                'sr':  np.ones(len(events), dtype='float'),
                'wmcr': leading_mu.id_sf,
                'tmcr': leading_mu.id_sf,
                'wecr': leading_e.id_sf,
                'tecr': leading_e.id_sf,
                'qcdcr': np.ones(len(events), dtype='float'),
            }
           
            ###
            # Reconstruction weights for electrons
            ###
                                                  
            reco = {
                'sr': np.ones(len(events), dtype='float'),
                'wmcr': np.ones(len(events), dtype='float'),
                'tmcr': np.ones(len(events), dtype='float'),
                'wecr': leading_e.reco_sf,
                'tecr': leading_e.reco_sf,
                'qcdcr': np.ones(len(events), dtype='float'),
            }

            ###
            # Isolation weights for muons
            ###

            '''
            isolation = {
                'sr': np.ones(len(events), dtype='float'),
                'wmcr': leading_mu.iso_sf,
                'tmcr': leading_mu.iso_sf,
                'wecr': np.ones(len(events), dtype='float'),
                'tecr': np.ones(len(events), dtype='float'),
                'qcdcr': np.ones(len(events), dtype='float'),
            }
            '''

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
            btagSFlight_uncorrelatedDown  = get_btag_weight('deepflav',self._year,'loose').btag_weight(
                j_iso.pt,
                j_iso.eta,
                j_iso.hadronFlavour,
                j_iso.isdflvL
            )

            if hasattr(events, "L1PreFiringWeight"): 
                weights.add('prefiring', events.L1PreFiringWeight.Nom, events.L1PreFiringWeight.Up, events.L1PreFiringWeight.Dn)
            weights.add('genw',events.genWeight)
            weights.add('nlo_qcd',nlo_qcd)
            weights.add('nlo_ewk',nlo_ewk)
            if 'cen' in nnlo_nlo:
                #weights.add('nnlo_nlo',nnlo_nlo['cen'])
                weights.add('qcd1',np.ones(len(events), dtype='float'), nnlo_nlo['qcd1up']/nnlo_nlo['cen'], nnlo_nlo['qcd1do']/nnlo_nlo['cen'])
                weights.add('qcd2',np.ones(len(events), dtype='float'), nnlo_nlo['qcd2up']/nnlo_nlo['cen'], nnlo_nlo['qcd2do']/nnlo_nlo['cen'])
                weights.add('qcd3',np.ones(len(events), dtype='float'), nnlo_nlo['qcd3up']/nnlo_nlo['cen'], nnlo_nlo['qcd3do']/nnlo_nlo['cen'])
                weights.add('ew1',np.ones(len(events), dtype='float'), nnlo_nlo['ew1up']/nnlo_nlo['cen'], nnlo_nlo['ew1do']/nnlo_nlo['cen'])
                weights.add('ew2G',np.ones(len(events), dtype='float'), nnlo_nlo['ew2Gup']/nnlo_nlo['cen'], nnlo_nlo['ew2Gdo']/nnlo_nlo['cen'])
                weights.add('ew3G',np.ones(len(events), dtype='float'), nnlo_nlo['ew3Gup']/nnlo_nlo['cen'], nnlo_nlo['ew3Gdo']/nnlo_nlo['cen'])
                weights.add('ew2W',np.ones(len(events), dtype='float'), nnlo_nlo['ew2Wup']/nnlo_nlo['cen'], nnlo_nlo['ew2Wdo']/nnlo_nlo['cen'])
                weights.add('ew3W',np.ones(len(events), dtype='float'), nnlo_nlo['ew3Wup']/nnlo_nlo['cen'], nnlo_nlo['ew3Wdo']/nnlo_nlo['cen'])
                weights.add('ew2Z',np.ones(len(events), dtype='float'), nnlo_nlo['ew2Zup']/nnlo_nlo['cen'], nnlo_nlo['ew2Zdo']/nnlo_nlo['cen'])
                weights.add('ew3Z',np.ones(len(events), dtype='float'), nnlo_nlo['ew3Zup']/nnlo_nlo['cen'], nnlo_nlo['ew3Zdo']/nnlo_nlo['cen'])
                weights.add('mix',np.ones(len(events), dtype='float'), nnlo_nlo['mixup']/nnlo_nlo['cen'], nnlo_nlo['mixdo']/nnlo_nlo['cen'])
                weights.add('muF',np.ones(len(events), dtype='float'), nnlo_nlo['muFup']/nnlo_nlo['cen'], nnlo_nlo['muFdo']/nnlo_nlo['cen'])
                weights.add('muR',np.ones(len(events), dtype='float'), nnlo_nlo['muRup']/nnlo_nlo['cen'], nnlo_nlo['muRdo']/nnlo_nlo['cen'])
            weights.add('pileup',pu)
            weights.add('trig', trig[region])
            weights.add('ids', ids[region])
            weights.add('reco', reco[region])
            #weights.add('isolation', isolation[region])
            weights.add('btagSF',btagSF)
            weights.add('btagSFbc_correlated',np.ones(len(events), dtype='float'), btagSFbc_correlatedUp/btagSF, btagSFbc_correlatedDown/btagSF)
            weights.add('btagSFbc_uncorrelated',np.ones(len(events), dtype='float'), btagSFbc_uncorrelatedUp/btagSF, btagSFbc_uncorrelatedDown/btagSF)
            weights.add('btagSFlight_correlated',np.ones(len(events), dtype='float'), btagSFlight_correlatedUp/btagSF, btagSFlight_correlatedDown/btagSF)
            weights.add('btagSFlight_uncorrelated',np.ones(len(events), dtype='float'), btagSFlight_uncorrelatedUp/btagSF, btagSFlight_uncorrelatedDown/btagSF)
            

        
        ###
        # Selections
        ###

        lumimask = np.ones(len(events), dtype='bool')
        if isData:
            lumimask = lumiMasks[self._year](events.run, events.luminosityBlock)
        selection.add('lumimask', lumimask)

        met_filters =  np.ones(len(events), dtype='bool')
        #if isData: met_filters = met_filters & events.Flag['eeBadScFilter']#this filter is recommended for data only
        for flag in AnalysisProcessor.met_filters[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters',met_filters)

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._met_triggers[self._year]:
            if not hasattr(events.HLT, path): continue
            triggers = triggers | events.HLT[path]
        selection.add('met_triggers', triggers)

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._singleelectron_triggers[self._year]:
            if not hasattr(events.HLT, path): continue
            triggers = triggers | events.HLT[path]
        selection.add('singleelectron_triggers', triggers)

        noHEMj = np.ones(len(events), dtype='bool')
        if self._year=='2018': noHEMj = (j_nHEM==0)

        noHEMmet = np.ones(len(events), dtype='bool')
        if self._year=='2018': noHEMmet = (met.pt>470)|(met.phi>-0.62)|(met.phi<-1.62)

        selection.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('isoneE', (e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
        #selection.add('leading_e_pt',(leading_e.pt>40))
        selection.add('noextrab', (j_ndflvL==0))
        selection.add('extrab', (j_ndflvL>0))
        selection.add('fatjet', (fj_nclean>0))
        selection.add('noHEMj', noHEMj)
        selection.add('noHEMmet', noHEMmet)
        selection.add('met100',(met.pt>100))
        selection.add('msd40',(leading_fj.msd_corr>40))
        selection.add('recoil_qcdcr', (u['qcdcr'].r>250))
        selection.add('mindphi_qcdcr', (ak.min(abs(u['qcdcr'].delta_phi(j_clean.T)), axis=1, mask_identity=False) < 0.1))
        selection.add('minDphi_qcdcr', (ak.min(abs(u['qcdcr'].delta_phi(fj_clean.T)), axis=1, mask_identity=False) > 1.5))
        selection.add('calo_qcdcr', ((abs(calomet.pt - met.pt) / u['qcdcr'].r)<0.5))

        regions = {
            'sr': ['msd40','fatjet', 'noHEMj','iszeroL','noextrab','met_filters','met_triggers','noHEMmet'],
            'wmcr': ['msd40','isoneM','fatjet','noextrab','noHEMj','met_filters','met_triggers'],
            'tmcr': ['msd40','isoneM','fatjet','extrab','noHEMj','met_filters','met_triggers'],
            'wecr': ['msd40','isoneE','fatjet','noextrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'tecr': ['msd40','isoneE','fatjet','extrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'qcdcr': ['recoil_qcdcr','mindphi_qcdcr','minDphi_qcdcr','calo_qcdcr','msd40','fatjet', 'noHEMj','iszeroL','noextrab','met_filters','met_triggers','noHEMmet'],
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
                  recoil=normalize(u[region].r, cut),
                  fjmass=normalize(leading_fj.msd_corr, cut),
                  TvsQCD=normalize(leading_fj.TvsQCD, cut),
                  weight=weight
            )
            if systematic is None:
                variables = {
                    'mindphirecoil':          ak.min(abs(u[region].delta_phi(j_clean.T)), axis=1,mask_identity=False),
                    'minDphirecoil':          ak.min(abs(u[region].delta_phi(fj_clean.T)), axis=1,mask_identity=False),
                    'CaloMinusPfOverRecoil':  abs(calomet.pt - met.pt) / u[region].r,
                    'met':                    met.pt,
                    'metphi':                 met.phi,
                    'mindphimet':             ak.min(abs(met.delta_phi(j_clean.T)), axis=1,mask_identity=False),
                    'minDphimet':             ak.min(abs(met.delta_phi(fj_clean.T)), axis=1,mask_identity=False),
                    'j1pt':                   leading_j.pt,
                    'j1eta':                  leading_j.eta,
                    'j1phi':                  leading_j.phi,
                    'fj1pt':                  leading_fj.pt,
                    'fj1eta':                 leading_fj.eta,
                    'fj1phi':                 leading_fj.phi,
                    'njets':                  j_nclean,
                    'ndflvL':                 j_ndflvL,
                    'nfjclean':               fj_nclean,
                }
                if region in mT:
                    variables['mT']           = mT[region]
                if 'e' in region:
                    variables['l1pt']      = leading_e.pt
                    variables['l1phi']     = leading_e.phi
                    variables['l1eta']     = leading_e.eta
                if 'm' in region:
                    variables['l1pt']      = leading_mu.pt
                    variables['l1phi']     = leading_mu.phi
                    variables['l1eta']     = leading_mu.eta
                for variable in output:
                    if variable not in variables:
                        continue
                    normalized_variable = {variable: normalize(variables[variable],cut)}
                    output[variable].fill(
                        region=region,
                        TvsQCD=normalize(leading_fj.TvsQCD,cut),
                        **normalized_variable,
                        weight=weight,
                    )
                output['TvsQCD'].fill(
                      region=region,
                      TvsQCD=normalize(leading_fj.TvsQCD, cut),
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

            if 'qcd' not in region:
                selection.add('recoil_'+region, (u[region].r>250))
                selection.add('mindphi_'+region, (ak.min(abs(u[region].delta_phi(j_clean.T)), axis=1, mask_identity=False) > 0.5))
                selection.add('minDphi_'+region, (ak.min(abs(u[region].delta_phi(fj_clean.T)), axis=1, mask_identity=False) > 1.5))
                selection.add('calo_'+region, ( (abs(calomet.pt - met.pt) / u[region].r) < 0.5))
                regions[region].insert(0, 'recoil_'+region)
                regions[region].insert(3, 'mindphi_'+region)
                regions[region].insert(4, 'minDphi_'+region)
                regions[region].insert(5, 'calo_'+region)
                if region in mT:
                    selection.add('mT_'+region, (mT[region]<150))

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
            output[key] *= scale
                
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
    ids         = load('data/ids.coffea')
    common      = load('data/common.coffea')

    processor_instance=AnalysisProcessor(year=options.year,
                                         xsec=xsec,
                                         corrections=corrections,
                                         ids=ids,
                                         common=common)

    save(processor_instance, 'data/hadmonotop'+options.name+'.processor')
