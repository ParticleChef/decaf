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
        #Values from https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2                                                      
        '2016postVFP': 36.31,
        '2016preVFP': 36.31,
        '2017': 41.48,
        '2018': 59.83
    }

    lumiMasks = {
        '2016postVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2016preVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2017': LumiMask("data/jsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
        '2018': LumiMask("data/jsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
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
        self._systematics = False
        self._skipJER = False

        self._samples = {
            'sr':('G1Jet','Z1Jets','Z2Jets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','SingleElectron','MET',
                    'TPhiTo2Chi_MPhi200_MChi50'),
            'wmcr':('WJets','Z1Jets','Z2Jets', 'G1Jet','DY','TT','ST','WW','WZ','ZZ','QCD','MET',
                    'TPhiTo2Chi_MPhi200_MChi50'),
            'tmcr':('WJets','Z1Jets','Z2Jets', 'G1Jet','DY','TT','ST','WW','WZ','ZZ','QCD','MET',
                    'TPhiTo2Chi_MPhi200_MChi50'),
            'wecr':('WJets','Z1Jets','Z2Jets', 'G1Jet','DY','TT','ST','WW','WZ','ZZ','QCD','SingleElectron','EGamma',
                    'TPhiTo2Chi_MPhi200_MChi50'),
            'tecr':('WJets','Z1Jets','Z2Jets', 'G1Jet','DY','TT','ST','WW','WZ','ZZ','QCD','SingleElectron','EGamma',
                    'TPhiTo2Chi_MPhi200_MChi50'),
            'zmcr':('WJets','Z1Jets','Z2Jets', 'G1Jet','DY','TT','ST','WW','WZ','ZZ','QCD','MET',
                    'TPhiTo2Chi_MPhi200_MChi50'),
            'zecr':('WJets','Z1Jets','Z2Jets', 'G1Jet','DY','TT','ST','WW','WZ','ZZ','QCD','SingleElectron','EGamma',
                    'TPhiTo2Chi_MPhi200_MChi50'),
            'gcr' :('WJets','Z1Jets','Z2Jets', 'G1Jet','DY','TT','ST','WW','WZ','ZZ','QCD','SinglePhoton','EGamma',
                    'TPhiTo2Chi_MPhi200_MChi50')
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
        self._singlephoton_triggers = {
            '2016': [
                'Photon175',
                'Photon165_HE10'
            ],
            '2017': [
                'Photon200'
            ],
            '2018': [
                'Photon200'
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
                'Photon200'
            ],
            '2018': [
                'Ele32_WPTight_Gsf',
                'Photon200'
            ]
        }
        self._singlemuon_triggers = {
            '2016': [
                'IsoMu24',
                'IsoTkMu24',
            ],
            '2017':
                [
                'IsoMu27',
            ],
            '2018':
                [
                'IsoMu24',
            ]
        }
        self._corrections = corrections
        self._ids = ids
        self._common = common

        #ptbins=np.arange(250,1210,38)

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
            'ut': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(65,350,1000, name='ut', label='U_{T}'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'uphi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='uphi', label='U_{phi}'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'met': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,600, name='met', label='MET $p_{T}$ [GeV]'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'metphi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='metphi', label='MET $\phi$'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'j1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(38,250,1200, name='j1pt', label='AK4 Leading Jet $p_{T}$'),
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
                hist.axis.Regular(38,250,1200, name='fj1pt', label='AK15 Leading SoftDrop Jet Pt'),
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
                hist.axis.Regular(38,250,1200, name='l1pt', label='Leading Lepton Pt'),
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
        fatjet_factory           = self._corrections['fatjet_factory']
        subjet_factory           = self._corrections['subjet_factory']
        met_factory              = self._corrections['met_factory']

        
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
        fatjets = fatjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiJet, events.fixedGridRhoFastjetAll), jec_cache)
        subjets = subjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiSubJet, events.fixedGridRhoFastjetAll), jec_cache)
        met = met_factory.build(events.MET, jets, {})

        shifts = [({"Jet": jets, "AK15PFPuppiSubJet": subjets, "AK15PFPuppiJet": fatjets, "MET": met}, None)]
        if self._systematics:
            shifts.extend([
                ({"Jet": jets.JES_jes.up, "AK15PFPuppiSubJet": subjets.JES_jes.up, "AK15PFPuppiSubJet": fatjets.JES_jes.up, "MET": met.JES_jes.up}, "JESUp"),
                ({"Jet": jets.JES_jes.down, "AK15PFPuppiSubJet": subjets.JES_jes.down, "AK15PFPuppiSubJet": fatjets.JES_jes.down, "MET": met.JES_jes.down}, "JESDown"),
                ({"Jet": jets, "AK15PFPuppiSubJet": subjets, "AK15PFPuppiSubJet": fatjets, "MET": met.MET_UnclusteredEnergy.up}, "UESUp"),
                ({"Jet": jets, "AK15PFPuppiSubJet": subjets, "AK15PFPuppiSubJet": fatjets, "MET": met.MET_UnclusteredEnergy.down}, "UESDown"),
            ])
            if not self._skipJER:
                shifts.extend([
                    ({"Jet": jets.JER.up, "AK15PFPuppiSubJet": subjets.JER.up, "AK15PFPuppiSubJet": fatjets.JER.up, "MET": met.JER.up}, "JERUp"),
                    ({"Jet": jets.JER.down, "AK15PFPuppiSubJet": subjets.JER.down, "AK15PFPuppiSubJet": fatjets.JER.down, "MET": met.JER.down}, "JERDown"),
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
        #get_mu_rochester_sf      = self._corrections['get_mu_rochester_sf'][self._year]
        get_met_xy_correction    = self._corrections['get_met_xy_correction']
        get_pu_weight            = self._corrections['get_pu_weight']    
        get_nlo_ewk_weight       = self._corrections['get_nlo_ewk_weight']    
        get_nnlo_nlo_weight      = self._corrections['get_nnlo_nlo_weight'][self._year]
        get_msd_corr             = self._corrections['get_msd_corr']
        get_btag_weight      = self._corrections['get_btag_weight']
        get_ttbar_weight     = self._corrections['get_ttbar_weight']
        
        isLooseElectron = self._ids['isLooseElectron'] 
        isTightElectron = self._ids['isTightElectron'] 
        isLooseMuon     = self._ids['isLooseMuon']     
        isTightMuon     = self._ids['isTightMuon']     
        #isLooseTau      = self._ids['isLooseTau']      
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
        # define leading mu
        leading_mu = ak.firsts(mu_tight)
        # define second mu for Z->mumu
        second_mu = ak.pad_none(mu_tight, target=2)[:,1]
        #dimu = ak.cartesian({"mu1": leading_mu, "mu2": second_mu}, nested=True)
        dimu = leading_mu + second_mu
        #dimu_mass = (dimu.mu1 + dimu.mu2).mass
        dimu_mass = dimu.mass

        e = events.Electron
## ReMOVE e cleanning from mu (eta-phi 0.3) 2024.4.1 (Not joke!)
        e['isclean'] = ak.all(e.metric_table(mu_loose) > -99999, axis=2)
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
        # define leading e
        leading_e = ak.firsts(e_tight)
        # define second e for Z->ee
        second_e = ak.pad_none(e_tight, target=2)[:,1]
        #diele = ak.cartesian({"e1": leading_e, "e2": second_e}, nested=True)
        diele = leading_e + second_e
        #diele_mass = (diele.e1 + diele.e2).mass
        diele_mass = diele.mass

        pho = events.Photon
        pho['isclean']=(
            ak.all(pho.metric_table(mu_loose) > 0.5, axis=2)
            & ak.all(pho.metric_table(e_loose) > 0.5, axis=2)
        )
        pho['T'] = ak.zip(
            {
                "r": pho.pt,
                "phi": pho.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        pho['isloose']=isLoosePhoton(pho,self._year)
        pho['istight']=isTightPhoton(pho,self._year)
        pho_clean=pho[pho.isclean]
        pho_loose=pho_clean[pho_clean.isloose]
        pho_tight=pho_clean[pho_clean.istight]
        pho_ntot=ak.num(pho, axis=1)
        pho_nloose=ak.num(pho_loose, axis=1)
        pho_ntight=ak.num(pho_tight, axis=1)
        leading_pho = ak.firsts(pho_tight)

        fj = events.AK15PFPuppiJet
        fj['pt'] = fj.subjets.sum().pt
        fj['msd_corr'] = get_msd_corr(fj)
        fj['isclean'] = (
            ak.all(fj.metric_table(mu_loose) > 1.5, axis=2)
            & ak.all(fj.metric_table(e_loose) > 1.5, axis=2)
#            & ak.all(fj.metric_table(pho_loose) > 1.5, axis=2)
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
#            & ak.all(j.metric_table(pho_loose) > 0.4, axis=2)
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
            'zecr'  : met+leading_e.T+second_e.T,
            'zmcr'  : met+leading_mu.T+second_mu.T,
            'gcr'   : met+leading_pho.T,
        }

        mT = {
            'sr'    : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            'wecr'  : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            'tecr'  : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            'wmcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
            'tmcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
            #'zecr'  : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            #'zmcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
            #'gcr'   : np.sqrt(2*leading_pho.pt*met.pt*(1-np.cos(met.delta_phi(leading_pho.T)))),
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
            if('TTTo' in dataset): 
                nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt) * get_ttbar_weight(genTops[:,1].pt))
                
            gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isA'] = (abs(gen.pdgId)==22) & gen.hasFlags(['isPrompt', 'fromHardProcess', 'isLastCopy']) & (gen.status == 1)
            
            genWs = gen[gen.isW] 
            genZs = gen[gen.isZ]
            genDYs = gen[gen.isZ&(gen.mass>30)]
            genAs = gen[gen.isA & (gen.pt > 100)]
            
            nlo_ewk = np.ones(len(events), dtype='float')
            if('WJets' in dataset): 
                #nlo_ewk = get_nlo_ewk_weight['w'](genWs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['w'](ak.max(genWs.pt, axis=1))
            elif('DY' in dataset): 
                nlo_ewk = get_nlo_ewk_weight['dy'](ak.max(genDYs.pt, axis=1))
            elif('Z1Jets' in dataset or 'Z2Jets' in dataset): 
                nlo_ewk = get_nlo_ewk_weight['z'](ak.max(genZs.pt, axis=1))
            elif('G1Jet' in dataset):
                nlo_ewk = get_nlo_ewk_weight['a'](ak.max(genAs.pt, axis=1))

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
                'zmcr': get_met_trig_weight(self._year, met.pt),
                'zecr': get_ele_trig_weight(self._year, leading_e.eta+leading_e.deltaEtaSC, leading_e.pt),
                'gcr': np.ones(len(events), dtype='float'),
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
                'zmcr': np.ones(len(events), dtype='float'),
                'zecr': np.ones(len(events), dtype='float'),
                'gcr':  np.ones(len(events), dtype='float'),
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
                'zmcr': np.ones(len(events), dtype='float'),
                'zecr': leading_e.reco_sf,
                'gcr':  np.ones(len(events), dtype='float'),
            }

            ###
            # Isolation weights for muons
            ###
            isolation = {
                'sr': np.ones(len(events), dtype='float'),
                'wmcr': leading_mu.iso_sf,
                'tmcr': leading_mu.iso_sf,
                'wecr': np.ones(len(events), dtype='float'),
                'tecr': np.ones(len(events), dtype='float'),
                'zmcr': leading_mu.iso_sf,
                'zecr': np.ones(len(events), dtype='float'),
                'gcr':  np.ones(len(events), dtype='float'),
            }

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
            weights.add('nlo_ewk',nlo_ewk)
            weights.add('pileup',pu)
            weights.add('trig', trig[region])
            weights.add('ids', ids[region])
            weights.add('reco', reco[region])
            weights.add('isolation', isolation[region])
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
            lumimask = AnalysisProcessor.lumiMasks[self._year](events.run, events.luminosityBlock)
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

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._singlephoton_triggers[self._year]:
            if not hasattr(events.HLT, path): continue
            triggers = triggers | events.HLT[path]
        selection.add('single_photon_triggers', ak.to_numpy(triggers))

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._singlemuon_triggers[self._year]:
            if path not in events.HLT.fields:
                continue
            triggers = triggers | events.HLT[path]
        selection.add('single_muon_triggers', ak.to_numpy(triggers))

        noHEMj = np.ones(len(events), dtype='bool')
        if self._year=='2018':
            noHEMj = (j_nHEM==0)

        noHEMmet = np.ones(len(events), dtype='bool')
        if self._year=='2018':
            noHEMmet = (met.pt>470)|(met.phi>-0.62)|(met.phi<-1.62)

        if ('WJetsToLNu' in dataset) & ('Pt' in dataset):
            remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt >120) ## W
            selection.add("exclude_wjets_greater_400", ak.to_numpy(ak.all(remove_overlap, axis=1)))
        else:
            selection.add("exclude_wjets_greater_400", np.full(len(events), True))

        if ('WJetsToLNu' in dataset) & (not ('Pt' in dataset)):
            remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt <= 120) ## w
            selection.add("exclude_wjets_less_400", ak.to_numpy(ak.all(remove_overlap, axis=1)))
        else:
            selection.add("exclude_wjets_less_400", np.full(len(events), True))

        selection.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(pho_nloose==0))
        selection.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(pho_nloose==0))
        selection.add('isoneE', (e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(pho_nloose==0))
        selection.add('isoneG', (e_nloose==0)&(mu_nloose==0)&(pho_nloose==1),(pho_ntight==1))
        selection.add('istwoM', (e_nloose==0)&(mu_ntight==2)&(mu_nloose==2)&(pho_nloose==0))
        selection.add('istwoE', (e_ntight==2)&(e_nloose==2)&(mu_nloose==0)&(pho_nloose==0))
        selection.add('one_ak4', (j_nclean>0))
        selection.add('one_ak15', (fj_nclean>0))
        selection.add('leading_fj250', (leading_fj.pt>250))

        selection.add('noextrab', (j_ndflvL==0))
        selection.add('extrab', (j_ndflvL>0))
        selection.add('oneb', (j_ndflvL==1))
        selection.add('noHEMj', noHEMj)
        selection.add('noHEMmet', noHEMmet)

        selection.add('met120',(met.pt<120))
        selection.add('met150',(met.pt>150))
        selection.add('diele60',(diele_mass>60))
        selection.add('diele120',(diele_mass<120))
        selection.add('dimu60',(dimu_mass>60))
        selection.add('dimu120',(dimu_mass<120))
        selection.add('leading_ele40',(leading_e.pt>40))

        regions = {
            'sr': [
                    'met_filters', 'met_triggers',
                    'noHEMj', 'noHEMmet',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'iszeroL',
                    'noextrab',
            ],
            'wmcr': [
                    'met_filters', 'met_triggers',
                    'noHEMj',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneM',
                    'noextrab',
                    'met150',
            ],
            'wecr': [
                    'met_filters', 'singleelectron_triggers',
                    'noHEMj',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneE',
                    'noextrab',
                    'met150',
            ],
            'tmcr': [
                    'met_filters', 'met_triggers',
                    'noHEMj',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneM',
                    'one_ak4',
                    'oneb',
                    'met150',
            ],
            'tecr': [
                    'met_filters', 'singleelectron_triggers',
                    'noHEMj',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneE',
                    'one_ak4',
                    'oneb',
                    'met150',
            ],
            'zmcr': [
                    'met_filters', 'met_triggers',
                    'noHEMj',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'istwoM',
                    'met120',
                    'dimu60', 'dimu120',
            ],
            'zecr': [
                    'met_filters', 'singleelectron_triggers',
                    'noHEMj',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'istwoE',
                    'met120',
                    'diele60', 'diele120',
                    'leading_ele40',
            ],
            'gcr': [
                    'met_filters', 'single_photon_triggers',
                    'noHEMj',
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
                    'ut' :                   u[region].r,
                    'uphi' :                 u[region].phi,
                    'met':                    met.pt,
                    'metphi':                 met.phi,
                    'j1pt':                   leading_j.pt,
                    'j1eta':                  leading_j.eta,
                    'j1phi':                  leading_j.phi,
                    'fj1pt':                  leading_fj.pt,
                    'fj1eta':                 leading_fj.eta,
                    'fj1phi':                 leading_fj.phi,
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
                if 'g' in region:
                    variables['l1pt']      = leading_pho.pt
                    variables['l1phi']     = leading_pho.phi
                    variables['l1eta']     = leading_pho.eta
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
                selection.add('recoil_'+region, (u[region].r>350))
                selection.add('mindphi_'+region, (ak.min(abs(u['sr'].delta_phi(j_clean.T)), axis=1, mask_identity=False) > 0.5))
                selection.add('minDphi_'+region, (ak.min(abs(u[region].delta_phi(fj_clean.T)), axis=1, mask_identity=False) > 1.5))
                #selection.add('calo_'+region, ( (abs(calomet.pt - met.pt) / u[region].r) < 0.5))
                regions[region].insert(0, 'recoil_'+region)
                regions[region].insert(3, 'mindphi_'+region)
                regions[region].insert(4, 'minDphi_'+region)
                #regions[region].insert(5, 'calo_'+region)
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
