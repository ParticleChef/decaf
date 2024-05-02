import numpy as np
from coffea.util import save
import awkward as ak

######
## Electron
## Electron_cutBased Int_t cut-based ID Fall17 V2
## (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
######


def isLooseElectron(e, year):
    
    if '2016' in year:
        year='2016'

    pt=e.pt
    eta=e.eta+e.deltaEtaSC
    dxy=e.dxy
    dz=e.dz
    loose_id=e.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2016":
        mask = (
            (pt > 10)
            & (abs(eta) < 1.4442)
#            & (abs(dxy) < 0.05)
#            & (abs(dz) < 0.1)
            & (loose_id >= 2)
        ) | (
            (pt > 10)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
#            & (abs(dxy) < 0.1)
#            & (abs(dz) < 0.2)
            & (loose_id >= 2)
        )
    elif year == "2017":
        mask = (
            (pt > 10)
            & (abs(eta) < 1.4442)
#            & (abs(dxy) < 0.05)
#            & (abs(dz) < 0.1)
            & (loose_id >= 2)
        ) | (
            (pt > 10)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
#            & (abs(dxy) < 0.1)
#            & (abs(dz) < 0.2)
            & (loose_id >= 2)
        )
    elif year == "2018":
        mask = (
            (pt > 10)
            & (abs(eta) < 1.4442)
#            & (abs(dxy) < 0.05)
#            & (abs(dz) < 0.1)
            & (loose_id >= 2)
        ) | (
            (pt > 10)
            & (abs(eta) > 1.5660)
            & (abs(eta) < 2.5)
#            & (abs(dxy) < 0.1)
#            & (abs(dz) < 0.2)
            & (loose_id >= 2)
        )
    return mask


def isTightElectron(e, year):
    if '2016' in year:
        year='2016'

    pt=e.pt
    eta=e.eta+e.deltaEtaSC
    dxy=e.dxy
    dz=e.dz
    tight_id=e.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2016":  # Trigger: HLT_Ele27_WPTight_Gsf_v
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
    elif year == "2017":  # Trigger: HLT_Ele35_WPTight_Gsf_v
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
    elif year == "2018":  # Trigger: HLT_Ele32_WPTight_Gsf_v
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
## Muon ID WPs:
## https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_selectors_Since_9_4_X
## Muon isolation WPs:
## https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection#Muon_Isolation
#######


def isLooseMuon(mu, year):
    
    if '2016' in year:
        year='2016'
        
    pt=mu.pt
    eta=mu.eta
    iso=mu.pfRelIso04_all
    loose_id=mu.looseId
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2016":
        mask = (pt > 20) & (abs(eta) < 2.4) & loose_id & (iso < 0.25)
    elif year == "2017":
        mask = (pt > 20) & (abs(eta) < 2.4) & loose_id & (iso < 0.25)
    elif year == "2018":
        mask = (pt > 15) & (abs(eta) < 2.4) & loose_id & (iso < 0.25)
    return mask


def isTightMuon(mu, year):
    
    if '2016' in year:
        year='2016'
        
    pt=mu.pt
    eta=mu.eta
    iso=mu.pfRelIso04_all
    tight_id=mu.tightId
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2016":
        mask = (pt > 30) & (abs(eta) < 2.4) & tight_id & (iso < 0.15)
    elif year == "2017":
        mask = (pt > 30) & (abs(eta) < 2.4) & tight_id & (iso < 0.15)
    elif year == "2018":
        mask = (pt > 30) & (abs(eta) < 2.4) & tight_id & (iso < 0.15)
    return mask


def isSoftMuon(mu, year):
    
    if '2016' in year:
        year='2016'
        
    pt=mu.pt
    eta=mu.eta
    iso=mu.pfRelIso04_all
    loose_id=mu.tightId
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2016":
        mask = (pt > 5) & (abs(eta) < 2.4) & tight_id & (iso > 0.15)
    elif year == "2017":
        mask = (pt > 5) & (abs(eta) < 2.4) & tight_id & (iso > 0.15)
    elif year == "2018":
        mask = (pt > 5) & (abs(eta) < 2.4) & tight_id & (iso > 0.15)
    return mask


######
## Tau
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
    
    if '2016' in year:
        year='2016'
        
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
    if year == "2016":
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
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
## Photon_cutBased Int_t cut-based ID bitmap, Fall17V2,
## (0:fail, 1:loose, 2:medium, 3:tight)
## Note: Photon IDs are integers, not bit masks
######


def isLoosePhoton(pho, year):
    
    if '2016' in year:
        year='2016'

    pt=pho.pt
    eta=pho.eta
    loose_id=pho.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2016":
        mask = (
            (pt > 20)
            & ~(abs(eta) > 1.4442)
            & (abs(eta) < 1.5660)
            & (abs(eta) < 2.5)
            & (loose_id >= 1)
        )
    elif year == "2017":
        mask = (
            (pt > 20)
            & ~(abs(eta) > 1.4442)
            & (abs(eta) < 1.5660)
            & (abs(eta) < 2.5)
            & (loose_id >= 1)
        )
    elif year == "2018":
        mask = (
            (pt > 20)
            & (~(abs(eta) > 1.4442) | (abs(eta) > 1.5660))
            & (abs(eta) < 2.5)
            & (loose_id >= 1)
        )
    return mask&(pho.electronVeto)


def isTightPhoton(pho, year):
    if '2016' in year:
        year='2016'

    pt=pho.pt
    eta=pho.eta
    tight_id=pho.cutBased
    
    mask = ~np.isnan(ak.ones_like(pt))
    if year == "2016":
        mask = (pt > 200) & (tight_id == 2) & ( abs(eta) < 1.479)
    elif year == "2017":
        mask = (pt > 230) & (tight_id == 2) & ( abs(eta) < 1.479)
    elif year == "2018":
        mask = (pt > 230) & (tight_id >= 2) & ( abs(eta) < 1.479)
    return mask&(pho.electronVeto)&(pho.isScEtaEB) #tight photons are barrel only


######
## Fatjet
## https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
## Tight working point including lepton veto (TightLepVeto)
######


def isGoodAK15(fj):
    
    pt=fj.pt
    eta=fj.eta
    jet_id=fj.jetId
    nhf=fj.neHEF
    chf=fj.chHEF
    
    mask = (
        (pt > 160) & (abs(eta) < 2.4) & ((jet_id & 6) == 6 ) & (nhf < 0.8) & (chf > 0.1)
    )
    return mask


######
## Jet
## https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
## Tight working point including lepton veto (TightLepVeto)
##
## For Jet ID flags, bit1 is Loose (always false in 2017 since it does not
## exist), bit2 is Tight, bit3 is TightLepVeto. The POG recommendation is to
## use Tight Jet ID as the standard Jet ID.
######
## PileupJetID
## https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
## Using Loose Pileup ID
##
## Note: There is a bug in 2016 UL in which bit values for Loose and Tight Jet
## Pileup IDs are accidentally flipped relative to 2017 UL and 2018 UL.
##
## For 2016 UL,
## Jet_puId = (passtightID*4 + passmediumID*2 + passlooseID*1).
##
## For 2017 UL and 2018 UL,
## Jet_puId = (passlooseID*4 + passmediumID*2 + passtightID*1).
######


def isGoodAK4(j, year):
    if '2016' in year:
        year='2016'
    
    pt=j.pt
    eta=j.eta
    jet_id=j.jetId
    pu_id=j.puId
    nhf=j.neHEF
    chf=j.chHEF
    
    mask = (pt > 30) & (abs(eta) < 2.4) & ((jet_id & 6) == 6)
    if year == "2016":
        mask = ((pt >= 50) & mask) | ((pt < 50) & mask & ((pu_id & 1) == 1))# & (nhf < 0.8) & (chf > 0.1)
    elif year == "2017":
        mask = ((pt >= 50) & mask) | ((pt < 50) & mask & ((pu_id & 4) == 4))# & (nhf < 0.8) & (chf > 0.1)
    elif year == "2018":
        mask = ((pt >= 50) & mask) | ((pt < 50) & mask & ((pu_id & 4) == 4))# & (nhf < 0.8) & (chf > 0.1)
    return mask


######
## HEM
######


def isHEMJet(j):

    pt=j.pt
    eta=j.eta
    phi=j.phi
    
    mask = (pt > 15) & (eta > -3.0) & (eta < -1.3) & (phi > -1.57) & (phi < -0.87)
    return mask


ids = {}
ids["isLooseElectron"] = isLooseElectron
ids["isTightElectron"] = isTightElectron
ids["isLooseMuon"] = isLooseMuon
ids["isTightMuon"] = isTightMuon
ids["isSoftMuon"] = isSoftMuon
ids["isLooseTau"] = isLooseTau
ids["isLoosePhoton"] = isLoosePhoton
ids["isTightPhoton"] = isTightPhoton
ids["isGoodAK4"] = isGoodAK4
ids["isGoodAK15"] = isGoodAK15
ids["isHEMJet"] = isHEMJet
save(ids, "data/ids.coffea")
