import numpy as np
import awkward as ak
from coffea.util import save
from coffea.nanoevents.methods import vector as v
from coffea import lookup_tools, jetmet_tools, util
from coffea.lookup_tools import extractor, dense_lookup
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory
import correctionlib
from correctionlib import convert

### Iso Track ID ###
def isTrackElectron(track, met_pt, met_phi, year):
    pt = track.pt
    eta = track.eta
    pdgId = abs(track.pdgId)
    relIso = track.pfRelIso03_all
    mT = np.sqrt(
        2 * pt * met_pt * (1 - np.cos(track.phi - met_phi))
    )
    mask = (
            (pt > 5)
            & (abs(eta) < 2.5)
            & (pdgId == 11)
            & (relIso < 0.2)
            & (mT < 100)
        )
    return mask

def isTrackMuon(track, met_pt, met_phi, year):
    pt = track.pt
    eta = track.eta
    pdgId = abs(track.pdgId)
    relIso = track.pfRelIso03_all
    mT = np.sqrt(
        2 * pt * met_pt * (1 - np.cos(track.phi - met_phi))
    )
    mask = (
            (pt > 5)
            & (abs(eta) < 2.5)
            & (pdgId == 13)
            & (relIso < 0.2)
            & (mT < 100)
        )
    return mask

def isTrackPion(track, met_pt, met_phi, year):
    pt = track.pt
    eta = track.eta
    pdgId = abs(track.pdgId)
    relIso = track.pfRelIso03_all
    mT = np.sqrt(
        2 * pt * met_pt * (1 - np.cos(track.phi - met_phi))
    )
    mask = (
            (pt > 10)
            & (abs(eta) < 2.5)
            & (pdgId == 211)
            & (relIso < 0.1)
            & (mT < 100)
        )
    return mask

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
            (pt > 220)
            & (abs(eta) < 1.4442)
            & (cutBased >= 3)
        ) | (
            (pt > 220)
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

def isGoodFatJet(jet, year):
    pt = jet.pt
    eta = jet.eta
    phi = jet.phi
    msd = jet.msoftdrop
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
        corr = evaluator["AK8PUPPI_TightLeptonVeto"]
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
    mask = (pt > 200) & (abs(eta) < 2.0) & (jetId == 1) & (vetoMap == 0) & (msd > 60)
    return mask

ids = {}
ids['isTrackElectron'] = isTrackElectron
ids['isTrackMuon'] = isTrackMuon
ids['isTrackPion'] = isTrackPion
ids['isVetoElectron'] = isVetoElectron
ids['isMediumElectron'] = isMediumElectron
ids['isLooseMuon'] = isLooseMuon
ids['isMediumMuon'] = isMediumMuon
ids['isMediumPhoton'] = isMediumPhoton
ids['isMediumTau'] = isMediumTau
ids['isGoodJet'] = isGoodJet
ids['isGoodFatJet'] = isGoodFatJet
save(ids, 'data/ids.coffea')
