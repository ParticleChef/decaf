import numpy as np
from coffea.util import save
import awkward as ak

######
## Electron
## cut-based ID RunIII Winter22
## (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun3
## https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIIRecommendations
######

def isLooseElectron(e, year):

    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'

    pt=e.pt
    eta=e.eta+e.deltaEtaSC
    #dxy=e.dxy ## (abs(dxy) < 0.05), (abs(dxy) < 0.1)
    #dz=e.dz ## (abs(dz) < 0.1), (abs(dz) < 0.2)
    loose_id=e.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":
        mask = (
            (pt > 10)
            & (abs(eta) < 1.4442)
            & (loose_id >= 1)
        ) | (
            (pt > 10)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
            & (loose_id >= 1)
        )
    elif year == "2023":
        mask = (
            (pt > 10)
            & (abs(eta) < 1.4442)
            & (loose_id >= 1)
        ) | (
            (pt > 10)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
            & (loose_id >= 1) # original is 2
        )
    return mask


def isTightElectron(e, year):
    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'

    pt=e.pt
    eta=e.eta+e.deltaEtaSC
    #dxy=e.dxy
    #dz=e.dz
    tight_id=e.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":  
        mask = (
            (pt > 40)
            & (abs(eta) < 1.4442)
#            & (abs(dxy) < 0.05)
#            & (abs(dz) < 0.1)
            & (tight_id == 4)
        ) | (
            (pt > 40)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.4)
#            & (abs(dxy) < 0.1)
#            & (abs(dz) < 0.2)
            & (tight_id == 4)
        )
    elif year == "2023":  
        mask = (
            (pt > 40)
            & (abs(eta) < 1.4442)
#            & (abs(dxy) < 0.05)
#            & (abs(dz) < 0.1)
            & (tight_id == 4)
        ) | (
            (pt > 40)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.4)
#            & (abs(dxy) < 0.1)
#            & (abs(dz) < 0.2)
            & (tight_id == 4)
        )
    return mask


#######
## Muon
## slimmedMuons after basic selection (pt > 15 || (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || 
##                                    passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt'))))
## Muon POG Recommendations:
## https://muon-wiki.docs.cern.ch/guidelines/recommendations/#medium-muon-id
## pfIsoId: (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)
## pfRelIso...
#######

def isLooseMuon(mu, year):
    
    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'
        
    pt=mu.pt
    eta=mu.eta
    #iso=mu.pfRelIso04_all
    #iso=mu.pfIsoId
    loose_id=mu.looseId
    isTracker=mu.isTracker
    ispfcan=mu.isPFcand
    isglobal=mu.isGlobal
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":
        mask = (pt > 10) & (abs(eta) < 2.4) & loose_id & isTracker & ispfcan & isglobal  #& (iso >= 2)
    else:
        mask = (pt > 10) & (abs(eta) < 2.4) & loose_id & isTracker & ispfcan & isglobal
    return mask


def isTightMuon(mu, year):
    
    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'
        
    pt=mu.pt
    eta=mu.eta
    #iso=mu.pfIsoId
    iso=mu.pfRelIso04_all # (iso < 0.15)
    tight_id=mu.tightId
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":
        mask = (pt > 30) & (abs(eta) < 2.4) & tight_id & (iso < 0.1)
    else:
        mask = (pt > 30) & (abs(eta) < 2.4) & tight_id & (iso < 0.1)
    return mask


def isSoftMuon(mu, year):
    
    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'
        
    pt=mu.pt
    eta=mu.eta
    iso=mu.pfRelIso04_all
    loose_id=mu.tightId
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":
        mask = (pt > 5) & (abs(eta) < 2.4) & tight_id & (iso > 0.15)
    else:
        mask = (pt > 5) & (abs(eta) < 2.4) & tight_id & (iso > 0.15)
    return mask


######
## Tau  ## NOT UPDATED FOR RUN3
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
## The decayModeFindingNewDMs: recommended for use with DeepTauv2p1, where decay
## modes 5 and 6 should be explicitly rejected.
## This should already be applied in NanoAOD.
##
## Tau_idDeepTau2017v2p1VSe ID working points (bitmask):
## 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose,
## 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
##
## Tau_idDeepTau2017v2p1VSjet ID working points (bitmask):
## 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose,
## 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
##
## Tau_idDeepTau2017v2p1VSmu ID working points (bitmask):
## 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight
######


def isLooseTau(tau, year):
    
    if '2022' in year:
        year='2022'
        
    pt = tau.pt
    eta = tau.eta
    ide = tau.idDeepTau2017v2p1VSe
    idj = tau.idDeepTau2017v2p1VSjet
    idm = tau.idDeepTau2017v2p1VSmu
    decayMode = tau.decayMode
    try:
        decayModeDMs=tau.decayModeFindingNewDMs
    except:
        decayModeDMs=~np.isnan(ak.ones_like(pt))

    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":
        mask = (
            (pt > 20)
            & (abs(eta) < 2.3)
            #& ~(decayMode == 5)
            #& ~(decayMode == 6)
            & decayModeDMs
            #& ((ide & 16) == 16)
            & ((idj & 4) == 4)
            #& ((idm & 2) == 2)
        )
    elif year == "2017":
        mask = (
            (pt > 20)
            & (abs(eta) < 2.3)
            #& ~(decayMode == 5)
            #& ~(decayMode == 6)
            & decayModeDMs
            #& ((ide & 16) == 16)
            & ((idj & 4) == 4)
            #& ((idm & 2) == 2)
        )
    elif year == "2018":
        mask = (
            (pt > 20)
            & (abs(eta) < 2.3)
            #& ~(decayMode == 5)
            #& ~(decayMode == 6)
            & decayModeDMs
            #& ((ide & 16) == 16)
            & ((idj & 4) == 4)
            #& ((idm & 2) == 2)
        )
    return mask


######
## Photon
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun3
## cut-based ID bitmap, RunIIIWinter22V1
## (0:fail, 1:loose, 2:medium, 3:tight)
## Note: Photon IDs are Unsigned integers, not bit masks
######


def isLoosePhoton(pho, year):
    
    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'

    pt=pho.pt
    eta=pho.eta
    loose_id=pho.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":
        mask = (
            (pt > 20)
            & (~(abs(eta) > 1.4442) | (abs(eta) > 1.5660))
            & (abs(eta) < 2.5)
            & (loose_id >= 1)
        )
    elif year == "2023":
        mask = (
            (pt > 20)
            & (~(abs(eta) > 1.4442) | (abs(eta) > 1.5660))
            & (abs(eta) < 2.5)
            & (loose_id >= 1)
        )
    return mask&(pho.electronVeto)


def isTightPhoton(pho, year):
    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'

    pt=pho.pt
    eta=pho.eta
    tight_id=pho.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2022":
        mask = (pt > 200) & (tight_id >= 2) & (~(abs(eta) > 1.4442) | (abs(eta) > 1.5660)) & (abs(eta) < 1.479)
    elif year == "2023":
        mask = (pt > 200) & (tight_id >= 2) & (~(abs(eta) > 1.4442) | (abs(eta) > 1.5660)) & (abs(eta) < 1.479)
    return mask&(pho.electronVeto)&(pho.isScEtaEB) #tight photons are barrel only


######
## Fatjet
## https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
## Tight working point including lepton veto (TightLepVeto)
######


def isGoodAK15(fj):
    
    pt=fj.pt
    eta=fj.eta
    jet_id=fj.jetId
    #nhf=fj.neHEF
    #chf=fj.chHEF
    
    mask = (
        (pt > 160) & (abs(eta) < 2.4) & ((jet_id & 6) == 6 )# & (nhf < 0.8) & (chf > 0.1)
    )
    return mask


######
## Jet
## slimmedJetsPuppi, i.e. ak4 PFJets Puppi with JECs applied, after basic selection (pt > 15)
##
## https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
## Jet ID flag: bit2 is tight, bit3 is tightLepVeto
##
## https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDRun3
## Default AK4 jets arre PUPPI jets.
## Analyses mainly using low pT (<50 GeV) jets or forward jets will profit the most from the PU JetID.
## Note that PileUpJetID should only be applied to JEC corrected jets with pT < 50 GeV and not applied above.
######


def isGoodAK4(j, year):
    if '2022' in year:
        year='2022'
    if '2023' in year:
        year='2023'
    
    pt=j.pt
    eta=j.eta
    jet_id=j.jetId
    #pu_id=j.puId
    nhf=j.neHEF
    chf=j.chHEF
    
    mask = (pt > 30) & (abs(eta) < 2.4) & ((jet_id & 6) == 6)
    if year == "2022":
        mask = ((pt >= 50) & mask) | ((pt < 50) & mask) # & ((pu_id & 1) == 1))# & (nhf < 0.8) & (chf > 0.1)
    else:
        mask = ((pt >= 50) & mask) | ((pt < 50) & mask) # & ((pu_id & 1) == 1))# & (nhf < 0.8) & (chf > 0.1)

    return mask



ids = {}
ids["isLooseElectron"] = isLooseElectron
ids["isTightElectron"] = isTightElectron
ids["isLooseMuon"] = isLooseMuon
ids["isTightMuon"] = isTightMuon
#ids["isSoftMuon"] = isSoftMuon
#ids["isLooseTau"] = isLooseTau
ids["isLoosePhoton"] = isLoosePhoton
ids["isTightPhoton"] = isTightPhoton
ids["isGoodAK4"] = isGoodAK4
ids["isGoodAK15"] = isGoodAK15
save(ids, "data/ids.coffea")
