# -*- coding: utf-8 -*-
import os, sys
import json
import uproot
import numpy as np
import awkward as ak
import tqdm
from multiprocessing import Pool, Manager
from coffea.nanoevents.methods import vector
#import vector
from coffea.util import load, save
import glob
import correctionlib
from correctionlib import convert

filepath = '/data/mc/RunIII2024Summer24NanoAODv15/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/*.root'
flist = glob.glob(filepath)[:100]
## Calcultate total file size
total_size = sum(os.path.getsize(f) for f in flist) / (1024**4)  # in TB
print(f'Total file size: {total_size:.2f} TB')

### Electron ID ###
def isVetoElectron(electron, year):
    pt = electron.pt
    eta = electron.eta #+ electron.deltaEtaSC
    cutBased = electron.cutBased
    miniIso = electron.miniPFRelIso_all
    mask = (
            (pt > 5)
            & (abs(eta) < 1.4442)
            & (cutBased >= 1)
            & (miniIso < 0.1)
        ) | (
            (pt > 5)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
            & (cutBased >= 1)
            & (miniIso < 0.1)
        )
    return mask

def isMediumElectron(electron, year):
    pt = electron.pt
    eta = electron.eta #+ electron.deltaEtaSC
    cutBased = electron.cutBased
    miniIso = electron.miniPFRelIso_all
    mask = (
            (pt > 10)
            & (abs(eta) < 1.4442)
            & (cutBased >= 3)
            & (miniIso < 0.1)
        ) | (
            (pt > 10)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
            & (cutBased >= 3)
            & (miniIso < 0.1)
        )
    return mask

### Muon ID ###
def isLooseMuon(muon, year):
    pt = muon.pt
    eta = muon.eta
    miniIso = muon.miniPFRelIso_all
    looseId = muon.looseId
    mask = (
            (pt > 5)
            & (abs(eta) < 2.4)
            & (looseId)
            & (miniIso < 0.2)
        )
    return mask

def isMediumMuon(muon, year):
    pt = muon.pt
    eta = muon.eta
    miniIso = muon.miniPFRelIso_all
    mediumId = muon.mediumId
    mask = (
            (pt > 10)
            & (abs(eta) < 2.4)
            & (mediumId)
            & (miniIso < 0.2)
        )
    return mask

### Photon ID ###
def isMediumPhoton(photon, year):
    pt = photon.pt
    eta = photon.eta
    cutBased = photon.cutBased
    mask = (
            (pt > 10)
            & (abs(eta) < 1.4442)
            & (cutBased >= 3)
        ) | (
            (pt > 10)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
            & (cutBased >= 3)
        )
    return mask

### Tau ID ###
def isMediumTau(tau, met_pt, met_phi, year):
    pt = tau.pt
    eta = tau.eta
    phi = tau.phi
    dz = tau.dz
    decayMode = tau.decayMode
    idj = tau.idDeepTau2018v2p5VSjet
    mT = np.sqrt(
        2 * pt * met_pt * (1 - np.cos(phi - met_phi))
    )
    mask = (
            (pt > 20)
            & (abs(eta) < 2.5)
            & (abs(dz) < 0.2)
            & ~(decayMode == 5)
            & ~(decayMode == 6)
            & (idj >= 5)
            & (mT < 100)
        )
    return mask

### Jet ID ###
''' deprecated
def isGoodJet(jet):
    pt = jet.pt
    eta = jet.eta
    jetId = jet.jetId
    mask = (pt > 30) & (abs(eta) < 2.4) & ((jetId & 6) == 6)
    return mask
'''
## https://cms-analysis-corrections.docs.cern.ch/corrections_era/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/JME/2025-07-17/#jetidjsongz
def isGoodJet(jet, year):
    pt = jet.pt
    eta = jet.eta
    phi = jet.phi
    chHEF = jet.chHEF
    neHEF = jet.neHEF
    chEmEF = jet.chEmEF
    neEmEF = jet.neEmEF
    muEF = jet.muEF
    chMultiplicity = jet.chMultiplicity
    neMultiplicity = jet.neMultiplicity
    multiplicity = chMultiplicity + neMultiplicity
    def getJetVeto(eta,phi):
        evaluator = correctionlib.CorrectionSet.from_file('data/JMESF/'+year+'/jetvetomaps.json.gz')
        corr = evaluator["Summer24Prompt24_RunBCDEFGHI_V1"]
        counts = ak.num(eta)
        eta, phi = ak.flatten(eta), ak.flatten(phi)
        out = corr.evaluate('jetvetomap',eta,phi)
        return ak.unflatten(out, counts)
    def getJetID(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity):
        evaluator = correctionlib.CorrectionSet.from_file('data/JMESF/'+year+'/jetid.json.gz')
        corr = evaluator["AK4PUPPI_TightLeptonVeto"]
        counts = ak.num(eta)
        eta, chHEF, neHEF, chEmEF, neEmEF, muEF = ak.flatten(eta), ak.flatten(chHEF), ak.flatten(neHEF), ak.flatten(chEmEF), ak.flatten(neEmEF), ak.flatten(muEF)
        chMultiplicity, neMultiplicity, multiplicity = ak.flatten(chMultiplicity), ak.flatten(neMultiplicity), ak.flatten(multiplicity)
        args = (
            eta,
            chHEF, neHEF, chEmEF, neEmEF, muEF,
            chMultiplicity, neMultiplicity, multiplicity,
        )
        out = corr.evaluate(*args)
        return ak.unflatten(out, counts)
    jetId = getJetID(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
    vetoMap = getJetVeto(eta, phi)
    mask = (pt > 30) & (abs(eta) < 2.4) & (jetId == 1) & (vetoMap == 0)
    return mask

def process_file(args):
    file = args
    nEvents = 0
    histo = {}

    with uproot.open(file, basketcache="uncompressed") as f:
        outname = file.split('/')[-1].split('.')[0].split('_')[-1]
        samplename = file.split('/')[-2]
        tr = f['Events']
        electron = tr.arrays(['Electron_pt', 'Electron_eta', 'Electron_cutBased', 'Electron_miniPFRelIso_all'], how='zip')['Electron']
        muon = tr.arrays(['Muon_pt', 'Muon_eta', 'Muon_looseId', 'Muon_mediumId', 'Muon_miniPFRelIso_all'], how='zip')['Muon']
        tau = tr.arrays(['Tau_pt', 'Tau_eta', 'Tau_phi', 'Tau_dz', 'Tau_decayMode', 'Tau_idDeepTau2018v2p5VSjet'], how='zip')['Tau']
        jet = tr.arrays(['Jet_pt', 'Jet_eta', 'Jet_phi','Jet_mass', 'Jet_chHEF', 'Jet_neHEF', 'Jet_chEmEF', 'Jet_neEmEF', 'Jet_muEF', 'Jet_chMultiplicity', 'Jet_neMultiplicity','Jet_btagUParTAK4B'], how='zip')['Jet']
        met = tr.arrays(['PuppiMET_pt', 'PuppiMET_phi'], how='zip')
        calomet = tr.arrays(['CaloMET_pt', 'CaloMET_phi'], how='zip')
        nEvents = len(tr)
        year = '2024'
        flat_variables = tr.arrays(['HLT_PFMET120_PFMHT120_IDTight','HLT_PFMET130_PFMHT130_IDTight','HLT_PFMET140_PFMHT140_IDTight','HLT_PFMETNoMu120_PFMHTNoMu120_IDTight','HLT_PFMETNoMu130_PFMHTNoMu130_IDTight','HLT_PFMETNoMu140_PFMHTNoMu140_IDTight',
        'Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_BadPFMuonDzFilter','Flag_eeBadScFilter','Flag_ecalBadCalibFilter'], how='zip')
        e_loose = electron[isVetoElectron(electron, year)]
        e_medium = electron[isMediumElectron(electron, year)]
        m_loose = muon[isLooseMuon(muon, year)]
        m_medium = muon[isMediumMuon(muon, year)]
        t_medium = tau[isMediumTau(tau, met.PuppiMET_pt, met.PuppiMET_phi, year)]
        j_good = jet[isGoodJet(jet, year)]
        b = j_good[j_good.btagUParTAK4B > 0.1272]
        n_e_loose = ak.num(e_loose, axis=1)
        n_e_medium = ak.num(e_medium, axis=1)
        n_m_loose = ak.num(m_loose, axis=1)
        n_m_medium = ak.num(m_medium, axis=1)
        n_t_medium = ak.num(t_medium, axis=1)
        n_j_good = ak.num(j_good, axis=1)
        n_b = ak.num(b, axis=1)

        nEvents = len(met.PuppiMET_pt)
        selection = (
            (met.PuppiMET_pt > 250)
            & (
                (flat_variables.HLT_PFMET120_PFMHT120_IDTight)
                | (flat_variables.HLT_PFMET130_PFMHT130_IDTight)
                | (flat_variables.HLT_PFMET140_PFMHT140_IDTight)
                | (flat_variables.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight)
                | (flat_variables.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight)
                | (flat_variables.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight)
            )
            & (flat_variables.Flag_goodVertices)
            & (flat_variables.Flag_globalSuperTightHalo2016Filter)
            & (flat_variables.Flag_HBHENoiseFilter)
            & (flat_variables.Flag_HBHENoiseIsoFilter)
            & (flat_variables.Flag_EcalDeadCellTriggerPrimitiveFilter)
            & (flat_variables.Flag_BadPFMuonFilter)
            & (flat_variables.Flag_BadPFMuonDzFilter)
            & (flat_variables.Flag_eeBadScFilter)
            & (flat_variables.Flag_ecalBadCalibFilter)
            & (n_e_loose == 0)
            & (n_m_loose == 0)
            & (n_t_medium == 0)
            & (n_j_good >= 2)
            & (n_b >= 1)
        )
        electron, muon, tau, jet, met = electron[selection], muon[selection], tau[selection], jet[selection], met[selection]
        calomet = calomet[selection]
        j_good = j_good[selection]
        flat_variables = flat_variables[selection]
        if len(electron) == 0:
            return 0

        j1 = ak.firsts(j_good)
        j2 = ak.pad_none(j_good, target=2)[:,1]
        # scalar HT
        scalarHT = ak.sum(j_good.pt, axis=1)
        # delta phi between jets and met
        j1_vector = ak.zip({
            "pt": j1.pt,
            "eta": j1.eta,
            "phi": j1.phi,
            "mass": j1.mass
        },
        with_name = "PtEtaPhiMLorentzVector",
        behavior = vector.behavior,
        )
        j2_vector = ak.zip({
            "pt": j2.pt,
            "eta": j2.eta,
            "phi": j2.phi,
            "mass": j2.mass
        },
        with_name = "PtEtaPhiMLorentzVector",
        behavior = vector.behavior,
        )
        met_vector = ak.zip({
            "pt": met.PuppiMET_pt,
            "eta": 0,
            "phi": met.PuppiMET_phi,
            "mass": 0
        },
        with_name = "PtEtaPhiMLorentzVector",
        behavior = vector.behavior,
        )
        dphi_j1_met = abs(j1_vector.delta_phi(met_vector))
        dphi_j2_met = abs(j2_vector.delta_phi(met_vector))
        puppiovercalo = met.PuppiMET_pt / calomet.CaloMET_pt

        cut2 = (
            (dphi_j1_met > 0.5)
            & (dphi_j2_met > 1.5)
            & (scalarHT > 300)
            & (puppiovercalo < 5)
        )
        electron, muon, tau, jet, met = electron[cut2], muon[cut2], tau[cut2], jet[cut2], met[cut2]
        flat_variables = flat_variables[cut2]

        met_return = met.PuppiMET_pt
        ## saving met values after all cuts in text file
        fname = 'met_test.txt'
        with open(fname, 'a') as f:
            for val in met_return:
                f.write(f"{val}\n")


        return nEvents

def main(flist):
    with Pool(10) as pool:
        results = list(tqdm.tqdm(pool.imap(
            process_file, [(file) for file in flist]), total=len(flist)))
        # sum all the results
        total = sum(results)
        print("Total events:", total)

if __name__ == "__main__":
    main(flist)