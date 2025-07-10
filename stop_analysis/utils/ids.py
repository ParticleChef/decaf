import numpy as np
import awkward as ak
from coffea.util import save
from coffea.nanoevents.methods import vector as v

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
def isGoodJet(jet):
    pt = jet.pt
    eta = jet.eta
    jetId = jet.jetId
    mask = (pt > 30) & (abs(eta) < 2.4) & ((jetId & 6) == 6)
    return mask


ids = {}
ids['isVetoElectron'] = isVetoElectron
ids['isMediumElectron'] = isMediumElectron
ids['isLooseMuon'] = isLooseMuon
ids['isMediumMuon'] = isMediumMuon
ids['isMediumPhoton'] = isMediumPhoton
ids['isMediumTau'] = isMediumTau
ids['isGoodJet'] = isGoodJet
save(ids, 'data/ids.coffea')
