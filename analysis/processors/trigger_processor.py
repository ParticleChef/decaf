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

class AnalysisProcessor(processor.ProcessorABC):
    lumis = {
        '2022pre': 7.9804,
        '2022post': 26.6717,
        '2023pre': 17.794,
        '2023post': 9.451,
        '2024':  108.95,
        '2025': 0.0
    }
    lumiMasks = {
        '2022pre': LumiMask('data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json'),
        '2022post': LumiMask('data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json'),
        '2023pre': LumiMask('data/lumiMask/Cert_Collisions2023_366442_370790_Golden.json'),
        '2023post': LumiMask('data/lumiMask/Cert_Collisions2023_366442_370790_Golden.json'),
        '2024': LumiMask('data/lumiMask/Cert_Collisions2024_378981_386951_Golden.json'),
        '2025': None
    }
    metFilters = {
        '2022pre': [
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
        '2022post': [
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
        '2023pre': [
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
        '2023post': [
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
        '2024': [
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
        '2025': []
    }

    def __init__(self, year, xsec, corrections, ids, common):
        self._year = year
        self._lumi = 1000 * float(self.lumis[year])
        self._xsec = xsec

        self._corrections = corrections
        self._ids = ids
        self._common = common
        
        self._samples = {
            'cat1_signal_trigger': ('TTtoLNu2Q','EGamma'),
            'cat2_electron_trigger': ('TTto2L2Nu','JetMET'),
            'cat3_muon_trigger': ('TTto2L2Nu','JetMET'),
        }
        self._signal_triggers ={
            '2022pre': [
                'PFMET120_PFMHT120_IDTight',
                'PFMET130_PFMHT130_IDTight',
                'PFMET140_PFMHT140_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight',
                'PFMETNoMu130_PFMHTNoMu130_IDTight',
                'PFMETNoMu140_PFMHTNoMu140_IDTight'
            ],
            '2022post': [
                'PFMET120_PFMHT120_IDTight',
                'PFMET130_PFMHT130_IDTight',
                'PFMET140_PFMHT140_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight',
                'PFMETNoMu130_PFMHTNoMu130_IDTight',
                'PFMETNoMu140_PFMHTNoMu140_IDTight'
            ],
            '2023pre': [
                'PFMET120_PFMHT120_IDTight',
                'PFMET130_PFMHT130_IDTight',
                'PFMET140_PFMHT140_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight',
                'PFMETNoMu130_PFMHTNoMu130_IDTight',
                'PFMETNoMu140_PFMHTNoMu140_IDTight'
            ],
            '2023post': [
                'PFMET120_PFMHT120_IDTight',
                'PFMET130_PFMHT130_IDTight',
                'PFMET140_PFMHT140_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight',
                'PFMETNoMu130_PFMHTNoMu130_IDTight',
                'PFMETNoMu140_PFMHTNoMu140_IDTight'
            ],
            '2024': [
                'PFMET120_PFMHT120_IDTight',
                'PFMET130_PFMHT130_IDTight',
                'PFMET140_PFMHT140_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight',
                'PFMETNoMu130_PFMHTNoMu130_IDTight',
                'PFMETNoMu140_PFMHTNoMu140_IDTight'
            ],
            '2025': []
        }
        self._reference_triggers = {
            '2022pre': [
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf'
            ],
            '2022post': [
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf'
            ],
            '2023pre': [
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf'
            ],
            '2023post': [
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf'
            ],
            '2024': [
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf'
            ],
            '2025': []
        }
        self._electron_triggers = {
            '2022pre': [
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Ele135_CaloIdVT_GsfTrkIdT',
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'DoubleEle25_CaloIdL_MW',
                'DoubleEle27_CaloIdL_MW',
                'DoubleEle33_CaloIdL_MW'
            ],
            '2022post': [
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Ele135_CaloIdVT_GsfTrkIdT',
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'DoubleEle25_CaloIdL_MW',
                'DoubleEle27_CaloIdL_MW',
                'DoubleEle33_CaloIdL_MW'
            ],
            '2023pre': [
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Ele135_CaloIdVT_GsfTrkIdT',
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'DoubleEle25_CaloIdL_MW',
                'DoubleEle27_CaloIdL_MW',
                'DoubleEle33_CaloIdL_MW'
            ],
            '2023post': [
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Ele135_CaloIdVT_GsfTrkIdT',
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'DoubleEle25_CaloIdL_MW',
                'DoubleEle27_CaloIdL_MW',
                'DoubleEle33_CaloIdL_MW'
            ],
            '2024': [
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Ele135_CaloIdVT_GsfTrkIdT',
                'Ele30_WPTight_Gsf',
                'Ele32_WPTight_Gsf',
                'Ele35_WPTight_Gsf',
                'Ele38_WPTight_Gsf',
                'Ele40_WPTight_Gsf',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
                'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
                'DoubleEle25_CaloIdL_MW',
                'DoubleEle27_CaloIdL_MW',
                'DoubleEle33_CaloIdL_MW'
            ],
            '2025': []
        }
        self._muon_triggers = {
            '2022pre': [
                'IsoMu20',
                'IsoMu24',
                'IsoMu27',
                'IsoMu24_eta2p1',
                'Mu50',
                'Mu55'
            ],
            '2022post': [
                'IsoMu20',
                'IsoMu24',
                'IsoMu27',
                'IsoMu24_eta2p1',
                'Mu50',
                'Mu55'
            ],
            '2023pre': [
                'IsoMu20',
                'IsoMu24',
                'IsoMu27',
                'IsoMu24_eta2p1',
                'Mu50',
                'Mu55'
            ],
            '2023post': [
                'IsoMu20',
                'IsoMu24',
                'IsoMu27',
                'IsoMu24_eta2p1',
                'Mu50',
                'Mu55'
            ],
            '2024': [
                'IsoMu20',
                'IsoMu24',
                'IsoMu27',
                'IsoMu24_eta2p1',
                'Mu50',
                'Mu55'
            ],
            '2025': []
        }
        self._photon_triggers = {
            '2022pre': [
                'Photon175',
                'Photon200'
            ],
            '2022post': [
                'Photon175',
                'Photon200'
            ],
            '2023pre': [
                'Photon175',
                'Photon200'
            ],
            '2023post': [
                'Photon175',
                'Photon200'
            ],
            '2024': [
                'Photon175',
                'Photon200'
            ],
            '2025': []
        }
        self._ht_triggers = {
            '2022pre': [
                'PFHT180',
                'PFHT250',
                'PFHT350',
                'PFHT430',
                'PFHT510',
                'PFHT590',
                'PFHT680',
                'PFHT780',
                'PFHT890',
                'PFHT1050',
            ],
            '2022post': [
                'PFHT180',
                'PFHT250',
                'PFHT350',
                'PFHT430',
                'PFHT510',
                'PFHT590',
                'PFHT680',
                'PFHT780',
                'PFHT890',
                'PFHT1050',
            ],
            '2023pre': [
                'PFHT180',
                'PFHT250',
                'PFHT350',
                'PFHT430',
                'PFHT510',
                'PFHT590',
                'PFHT680',
                'PFHT780',
                'PFHT890',
                'PFHT1050',
            ],
            '2023post': [
                'PFHT180',
                'PFHT250',
                'PFHT350',
                'PFHT430',
                'PFHT510',
                'PFHT590',
                'PFHT680',
                'PFHT780',
                'PFHT890',
                'PFHT1050',
            ],
            '2024': [
                'PFHT180',
                'PFHT250',
                'PFHT350',
                'PFHT430',
                'PFHT510',
                'PFHT590',
                'PFHT680',
                'PFHT780',
                'PFHT890',
                'PFHT1050',
            ],
            '2025': []
        }

        self.make_output = lambda: {
            'sumw': 0.0,
            'template': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.StrCategory([], name='systematic', growth=True),
                hist.axis.Variable([250,260,270,280,290,300,350,400,500,800], name='met', label=r'$E_{T}^{miss}$ (GeV)'),
                storage=hist.storage.Weight(),
            ),
            'metpt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.StrCategory([], name='systematic', growth=True),
                hist.axis.IntCategory([0, 1], name='signal_trigger', label='Signal Trigger'),
                hist.axis.Variable([100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,350,400,500,800], name='met', label=r'$E_{T}^{miss}$ (GeV)'),
                storage=hist.storage.Weight(),
            ),
        }
    
    def process(self, events):
        isData = not hasattr(events, 'genWeight')
        selection = PackedSelection(dtype="uint64")
        weights = Weights(len(events), storeIndividual=True)
        output = self.make_output()
        if not isData:
            output['sumw'] = ak.sum(events.genWeight)
        
        dataset = events.metadata['dataset']
        selected_regions = []
        for region, samples in self._samples.items():
            for sample in samples:
                if sample in dataset:
                    selected_regions.append(region)
                    break
                    
        """ Initialize the physics objects """
        ### read the ids
        isTrackElectron = self._ids['isTrackElectron']
        isTrackMuon = self._ids['isTrackMuon']
        isTrackPion = self._ids['isTrackPion']
        isVetoElectron = self._ids['isVetoElectron']
        isMediumElectron = self._ids['isMediumElectron']
        isLooseMuon = self._ids['isLooseMuon']
        isMediumMuon = self._ids['isMediumMuon']
        isMediumPhoton = self._ids['isMediumPhoton']
        isMediumTau = self._ids['isMediumTau']
        isGoodJet = self._ids['isGoodJet']
        isGoodFatJet = self._ids['isGoodFatJet']

        ### read the corrections
        get_pu_weight = self._corrections['get_pu_weight']
        get_jec_correction = self._corrections['get_jec_correction']
        get_fjec_correction = self._corrections['get_fjec_correction']

        ### Initialize global quantities
        npv = events.PV.npvsGood
        run = events.run
        if '2022' in self._year or '2023' in self._year:
            met = events.MET
        else:
            met = events.PuppiMET
        calo_met = events.CaloMET

        ### Tracks
        trk = events.IsoTrack
        trk['isTrackElectron'] = isTrackElectron(trk, met.pt, met.phi, self._year)
        trk['isTrackMuon'] = isTrackMuon(trk, met.pt, met.phi, self._year)
        trk['isTrackPion'] = isTrackPion(trk, met.pt, met.phi, self._year)
        trk_e = trk[trk.isTrackElectron]
        trk_m = trk[trk.isTrackMuon]
        trk_pi = trk[trk.isTrackPion]

        ### Electrons
        e = events.Electron
        e['isveto'] = isVetoElectron(e, self._year)
        e['ismedium'] = isMediumElectron(e, self._year)
        e['T'] = ak.zip({
            'r': e.pt,
            'phi': e.phi,
        }, with_name='PolarTwoVector', behavior=vector.behavior)
        e_veto = e[e.isveto]
        e_medium = e[e.ismedium]
        #print(e_medium.pt)
        mT_e = np.sqrt(2 * e_veto.pt * met.pt * (1 - np.cos(met.delta_phi(e_veto.T))))
        
        ### Muons
        m = events.Muon
        m['isloose'] = isLooseMuon(m, self._year)
        m['ismedium'] = isMediumMuon(m, self._year)
        m['T'] = ak.zip({
            'r': m.pt,
            'phi': m.phi,
        }, with_name='PolarTwoVector', behavior=vector.behavior)
        m_loose = m[m.isloose]
        m_medium = m[m.ismedium]
        mT_m = np.sqrt(2 * m_loose.pt * met.pt * (1 - np.cos(met.delta_phi(m_loose.T))))

        ### Photons
        p = events.Photon
        p['ismedium'] = isMediumPhoton(p, self._year)
        p['T'] = ak.zip({
            'r': p.pt,
            'phi': p.phi,
        }, with_name='PolarTwoVector', behavior=vector.behavior)
        p_medium = p[p.ismedium]

        ### Taus
        t = events.Tau
        t['ismedium'] = isMediumTau(t, met.pt, met.phi, self._year)
        t['T'] = ak.zip({
            'r': t.pt,
            'phi': t.phi,
        }, with_name='PolarTwoVector', behavior=vector.behavior)
        t_medium = t[t.ismedium]
        #print(t_medium.pt)

        ### Jets
        j = events.Jet # Events/Jet_*
        ### Appling JECs
        rho_density = events.Rho.fixedGridRhoFastjetAll
        jec_corr = get_jec_correction(self._year, j.pt, j.eta, j.phi, rho_density, j.area, run, isData)
        j['pt'] = j.pt * jec_corr
        j['mass'] = j.mass * jec_corr

        ### Appling JetID
        j['isgood'] = isGoodJet(j, self._year)
        j['T'] = ak.zip({
            'r': j.pt,
            'phi': j.phi,
        }, with_name='PolarTwoVector', behavior=vector.behavior)
        ## b-tagging
        j['isupartL'] = (j.btagUParTAK4B>0.0246)
        j['isupartM'] = (j.btagUParTAK4B>0.1272)
        j['isupartT'] = (j.btagUParTAK4B>0.4648)

        j_good = j[j.isgood]
        j1 = ak.firsts(j_good)
        j2 = ak.pad_none(j_good, target=2)[:,1]
        b = j_good[j_good.isupartM]
        b1 = ak.firsts(b)
        b_tight = j_good[j_good.isupartT]
        b1_tight = ak.firsts(b_tight)
        b_loose = j_good[j_good.isupartL]
        b1_loose = ak.firsts(b_loose)
        nb_loose = ak.num(b_loose, axis=1)
        nb_medium = ak.num(b, axis=1)
        nb_tight = ak.num(b_tight, axis=1)


        ### FatJets
        fj = events.FatJet
        ### Appling JECs
        fjec_corr = get_fjec_correction(self._year, fj.pt, fj.eta, fj.phi, rho_density, fj.area, run, isData)

        fj['pt'] = fj.pt * fjec_corr
        fj['mass'] = fj.mass * fjec_corr
        ### Appling FatJetID

        fj['isgood'] = isGoodFatJet(fj, self._year)
        
        fj['T'] = ak.zip({
            'r': fj.pt,
            'phi': fj.phi,
        }, with_name='PolarTwoVector', behavior=vector.behavior)
        fj_good = fj[fj.isgood]
        nfj_good = ak.num(fj_good, axis=1)

        ### Scalar HT
        scalarHT = ak.sum(j_good.pt, axis=1)

        """ Variables for selection """
        ### lumimask
        lumimask = np.ones(len(events), dtype=bool)
        if isData:
            lumimask = AnalysisProcessor.lumiMasks[self._year](events.run, events.luminosityBlock)
        selection.add('lumimask', lumimask)

        ### met filters
        met_filters = np.ones(len(events), dtype=bool)
        for flag in AnalysisProcessor.metFilters[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters', met_filters)

        ### triggers
        signal_triggers = np.zeros(len(events), dtype=bool)
        for path in self._signal_triggers[self._year]:
            if not hasattr(events.HLT, path):
                continue
            signal_triggers = signal_triggers | events.HLT[path]
        selection.add('signal_trigger', signal_triggers)

        reference_triggers = np.zeros(len(events), dtype=bool)
        for path in self._reference_triggers[self._year]:
            if not hasattr(events.HLT, path):
                continue
            reference_triggers = reference_triggers | events.HLT[path]
        selection.add('reference_trigger', reference_triggers)

        ### Number of objects
        n_trk_e = ak.num(trk_e, axis=1)
        n_trk_m = ak.num(trk_m, axis=1)
        n_trk_pi = ak.num(trk_pi, axis=1)
        n_e_veto = ak.num(e_veto, axis=1)
        n_e_medium = ak.num(e_medium, axis=1)
        n_m_loose = ak.num(m_loose, axis=1)
        n_m_medium = ak.num(m_medium, axis=1)
        n_p_medium = ak.num(p_medium, axis=1)
        n_t_medium = ak.num(t_medium, axis=1)
        n_j_good = ak.num(j_good, axis=1)
        n_b = ak.num(b, axis=1)

        ### Opening angle between jets and MET
        j1_met_dphi = np.abs(j1.delta_phi(met))
        j2_met_dphi = np.abs(j2.delta_phi(met))

        ### Opening angle for third, fourth and fifth jets if present
        j3 = ak.pad_none(j_good, target=3)[:,2]
        j4 = ak.pad_none(j_good, target=4)[:,3]
        j5 = ak.pad_none(j_good, target=5)[:,4]

        j3_met_dphi = np.abs(j3.delta_phi(met))
        j4_met_dphi = np.abs(j4.delta_phi(met))
        j5_met_dphi = np.abs(j5.delta_phi(met))
        
        """ Define the selections """

        selection.add('zero_trk_e', n_trk_e == 0)
        selection.add('zero_trk_m', n_trk_m == 0)
        selection.add('zero_trk_pi', n_trk_pi == 0)
        selection.add('zero_e', n_e_veto == 0)
        selection.add('zero_m', n_m_loose == 0)
        selection.add('zero_t', n_t_medium == 0)
        selection.add('one_veto_lepton', ((n_e_veto == 1) & (n_m_loose == 0)) | ((n_e_veto == 0) & (n_m_loose == 1)))
        selection.add('one_e', n_e_medium >= 1)
        selection.add('one_m', n_m_medium == 1)
        selection.add('one_p', n_p_medium == 1)
        selection.add('zero_b', n_b == 0)
        selection.add('one_b', n_b >= 1)
        selection.add('one_b_tight', ak.num(b_tight, axis=1) >=1)
        selection.add('one_b_loose', ak.num(b_loose, axis=1) >=1)
        selection.add('exact_one_b', n_b == 1)
        selection.add('two_b', n_b >= 2)
        selection.add('two_j', n_j_good >= 2)
        selection.add('five_j', n_j_good >= 5)
        selection.add('ht_300', scalarHT > 300)
        selection.add('met_250', met.pt > 250)
        selection.add('mT_100', (ak.all(mT_m < 100, axis=1) & ak.all(mT_e < 100, axis=1)))
        selection.add('puppi/calo', met.pt / calo_met.pt < 5)
        selection.add('opening_angles_preselection', (j1_met_dphi > 0.5) & (j2_met_dphi > 0.15) & ak.fill_none(j3_met_dphi > 0.15, True))
        selection.add('opening_angles_highDeltaM', (j1_met_dphi > 0.5) & (j2_met_dphi > 0.5) & (j3_met_dphi > 0.5) & (j4_met_dphi > 0.5))
        selection.add('opening_angles_QCDCR', (j1_met_dphi < 0.5) | (j2_met_dphi < 0.15) | ak.fill_none(j3_met_dphi < 0.15, False))
        selection.add('opening_angles_QCDCR_highDeltaM', (j1_met_dphi < 0.5) | (j2_met_dphi < 0.5) | (j3_met_dphi < 0.5) | (j4_met_dphi < 0.5))

        regions = {
            'cat1_signal_trigger': [
                'lumimask', 'met_filters',
                'reference_trigger',
                'zero_trk_e', 'zero_trk_m', 'zero_trk_pi',
                'zero_m', 'zero_t', 'one_e', 'two_j',
                'puppi/calo',
                'ht_300', 'opening_angles_preselection'
            ],
            'cat2_electron_trigger': [
                'lumimask', 'met_filters',
                'signal_trigger',
                'zero_trk_e', 'zero_trk_m', 'zero_trk_pi',
                'zero_m', 'zero_t', 'one_e', 'two_j',
                'puppi/calo',
                'ht_300', 'opening_angles_preselection'
            ],
            'cat3_muon_trigger': [
                'lumimask', 'met_filters',
                'signal_trigger',
                'zero_trk_e', 'zero_trk_m', 'zero_trk_pi',
                'zero_e', 'zero_t', 'one_m', 'two_j',
                'puppi/calo',
                'ht_300', 'opening_angles_preselection'
            ],
        }
        if not isData:
            weights.add('genweight', events.genWeight)
            # PU reweighting
            pu_nom, pu_up, pu_down = get_pu_weight(self._year, events.Pileup.nTrueInt)
            weights.add('pileup', pu_nom)

        def normalize(val, cut):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar

        def fill(region, systematic):
            cut = selection.all(*regions[region])
            sys_name = 'nominal' if systematic is None else systematic
            if systematic in weights.variations:
                weight = weights.weight(modifier=systematic)[cut]
            else:
                weight = weights.weight()[cut]
            output['template'].fill(
                region=region,
                systematic=sys_name,
                met=met.pt[cut],
                weight=weight
            )
            output['metpt'].fill(
                region=region,
                systematic=sys_name,
                signal_trigger= signal_triggers[cut],
                met=met.pt[cut],
                weight=weight
            )
            if systematic is None:
                variables = {
                }
                for variable in output:
                    if variable not in variables:
                        continue
                    normalized_variable = {variable: normalize(variables[variable], cut)}
                    output[variable].fill(
                        region=region,
                        systematic=sys_name,
                        **normalized_variable,
                        weight=weight
                    )
        shift_name = None
        if shift_name is None:
            systematics = [None] + list(weights.variations)
        else:
            systematics = [shift_name]
        for region in regions:
            if region not in selected_regions:
                continue
            for systematic in systematics:
                if isData and systematic is not None:
                    continue
                fill(region, systematic)
        scale = 1
        if self._xsec[dataset] > 0:
            scale = self._lumi * self._xsec[dataset]
        
        for key in output:
            if key == 'sumw':
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

    save(processor_instance, 'data/stop_'+options.name+'.processor')