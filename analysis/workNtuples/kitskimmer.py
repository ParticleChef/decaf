import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt #type: ignore
import mplhep as hep
import glob
import os

moritz_global_path = '/data/mTopnTuples/ntuples_kit/'

### Load data
moritz_ntuples = glob.glob(moritz_global_path + '/MET_2018*/*.root')
print(moritz_ntuples)

original_events = 0
cut_events = 0
for idx, f in enumerate(moritz_ntuples):
    histo = {}
    f = uproot.open(f)
    tr = f["Events"]
    if idx == 0:
        print(tr.keys())
    Lumi = tr['Evt_Lumi'].array()
    Run = tr['Evt_Run'].array()
    Event = tr['Evt_ID'].array()
    PFMETNoMu120_PFHT60 = tr.arrays(["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"], how="zip")['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60']
    PFMETNoMu120 = tr.arrays(["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"], how="zip")['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight']
    single_photon_trigger = tr.arrays(["HLT_Photon200"], how="zip")['HLT_Photon200']
    Ele32_WPTight_Gsf = tr.arrays(["HLT_Ele32_WPTight_Gsf"], how="zip")['HLT_Ele32_WPTight_Gsf']
    LoosePhoton = tr.arrays(["LoosePhoton_Pt", "LoosePhoton_Eta", "LoosePhoton_Phi"], how="zip")['LoosePhoton']
    TightPhoton = tr.arrays(["TightPhoton_Pt", "TightPhoton_Eta", "TightPhoton_Phi"], how="zip")['TightPhoton']
    LooseElectron = tr.arrays(["LooseElectron_Pt", "LooseElectron_Eta","LooseElectron_EtaSC", "LooseElectron_Phi"], how="zip")['LooseElectron']
    TightElectron = tr.arrays(["TightElectron_Pt", "TightElectron_Eta","TightElectron_EtaSC", "TightElectron_Phi"], how="zip")['TightElectron']
    LooseMuon = tr.arrays(["LooseMuon_Pt", "LooseMuon_Eta", "LooseMuon_Phi"], how="zip")['LooseMuon']
    TightMuon = tr.arrays(["TightMuon_Pt", "TightMuon_Eta", "TightMuon_Phi"], how="zip")['TightMuon']
    metxy_pt = tr['MET_T1XY_Pt_nom'].array()
    metxy_phi = tr['MET_T1XY_Phi_nom'].array()
    met_pt = tr['MET_T1_Pt_nom'].array()
    met_phi = tr['MET_T1_Phi_nom'].array()
    genmet_pt = tr['GenMET_Pt_nom'].array()
    genmet_phi = tr['GenMET_Phi_nom'].array()
    puppimet_pt = tr['PuppiMET_Pt_nom'].array()
    puppimet_phi = tr['PuppiMET_Phi_nom'].array()
    Recoil = tr.arrays(["Hadr_Recoil_MET_T1XY_Pt_nom", "Hadr_Recoil_MET_T1XY_Phi_nom"], how="zip")
    Jets_outside_lead_AK15Jet_taggedL_nom = tr['Jets_outside_lead_AK15Jet_taggedL_nom'].array()
    ak15 = tr.arrays(["AK15Jet_Pt_nom", "AK15Jet_Eta_nom", "AK15Jet_Phi_nom", "AK15Jet_M_nom","AK15Jet_chHEF_nom","AK15Jet_neHEF_nom"], how="zip")['AK15Jet']
    nak15 = tr['N_AK15Jets_nom'].array()
    N_HEMJets_nom = tr['N_HEMJets_nom'].array()
    N_Jets_nom = tr['N_Jets_nom'].array()
    N_AK15Jets_nom = tr['N_AK15Jets_nom'].array()
    N_AK15SoftDropJets_nom = tr['N_AK15SoftDropJets_nom'].array()
    MT = tr['M_W_transverse_MET_T1XY_nom'].array()

    DeltaPhi_AK15Jet_Hadr_Recoil_MET_T1XY_nom = (tr['DeltaPhi_AK15Jet_Hadr_Recoil_MET_T1XY_nom'].array() > 1.5)
    ALL_DeltaPhi_AK4Jet_Hadr_Recoil_MET_T1XY_nom = (tr['DeltaPhi_AK4Jet_Hadr_Recoil_MET_T1XY_nom'].array())
    N_outofAK15 = ak.sum(Jets_outside_lead_AK15Jet_taggedL_nom, axis=1)
    met_trigger = (PFMETNoMu120_PFHT60 | PFMETNoMu120)
    single_electron_trigger = (single_photon_trigger | Ele32_WPTight_Gsf)

    
    npho_tight = ak.num(TightPhoton.Pt)
    npho_loose = ak.num(LoosePhoton.Pt)
    nmu_tight = ak.num(TightMuon.Pt)
    nmu_loose = ak.num(LooseMuon.Pt)
    nel_tight = ak.num(TightElectron.Pt)
    nel_loose = ak.num(LooseElectron.Pt)
    
    nhem_mask = (N_HEMJets_nom == 0)
    recoil_mask = (
        (Recoil.Hadr_Recoil_MET_T1XY_Pt_nom > 350)
        & ak.firsts(DeltaPhi_AK15Jet_Hadr_Recoil_MET_T1XY_nom)
        & (ak.min(ALL_DeltaPhi_AK4Jet_Hadr_Recoil_MET_T1XY_nom, axis=1) > 0.5)
    )
    ak15_mask = (
        (N_AK15Jets_nom >= 1)
        & (ak.firsts(ak15.Pt_nom) >= 250)
        & (N_Jets_nom >= 1)
        & (ak.firsts(ak15.chHEF_nom) > 0.1) & (ak.firsts(ak15.neHEF_nom) < 0.8)
        & (N_AK15SoftDropJets_nom == N_AK15Jets_nom)
    )
    gcr_mask = (
        (npho_loose == 1) & (npho_tight == 1)
        & (nmu_loose == 0)
        & (nel_loose == 0)
        & (N_outofAK15 == 0)
    )
    wecr_mask = (
        (npho_loose == 0) & (npho_tight == 0)
        & (nmu_loose == 0)
        & (nel_loose == 1) & (nel_tight == 1)
        & (N_outofAK15 == 0)
        & ak.fill_none(((ak.num(MT) == 1) & (ak.firsts(MT < 150))), False)
        & (met_pt > 150)
    )
    wmcr_mask = (
        (npho_loose == 0) & (npho_tight == 0)
        & (nel_loose == 0)
        & (nmu_loose == 1) & (nmu_tight == 1)
        & (N_outofAK15 == 0)
        & ak.fill_none(((ak.num(MT) == 1) & (ak.firsts(MT < 150))), False)
        & (met_pt > 150)
    )
    #mask = nhem_mask & ak15_mask & recoil_mask & single_photon_trigger & gcr_mask
    mask = nhem_mask & ak15_mask & recoil_mask & met_trigger & wmcr_mask

    # change None to False
    mask_final = ak.fill_none(mask, False)
    original_events += len(single_photon_trigger)
    cut_events += len(single_photon_trigger[mask_final])
    print(met_pt[mask_final])
    histo['LoosePhoton'] = LoosePhoton[mask_final]
    histo['TightPhoton'] = TightPhoton[mask_final]
    histo['LooseElectron'] = LooseElectron[mask_final]
    histo['TightElectron'] = TightElectron[mask_final]
    histo['LooseMuon'] = LooseMuon[mask_final]
    histo['TightMuon'] = TightMuon[mask_final]
    histo['FJet'] = ak15[mask_final]
    histo['Recoil'] = Recoil[mask_final]
    histo['MET_Pt'] = met_pt[mask_final]
    histo['MET_Phi'] = met_phi[mask_final]
    histo['METXY_Pt'] = metxy_pt[mask_final]
    histo['METXY_phi'] = metxy_phi[mask_final]
    histo['GenMET_Pt'] = genmet_pt[mask_final]
    histo['GenMET_Phi'] = genmet_phi[mask_final]
    histo['PuppiMET_Pt'] = puppimet_pt[mask_final]
    histo['PuppiMET_Phi'] = puppimet_phi[mask_final]
    histo['MT'] = MT[mask_final]
    histo['Lumi'] = Lumi[mask_final]
    histo['Run'] = Run[mask_final]
    histo['Event'] = Event[mask_final]
    with uproot.recreate(f"kit_wmcr/kit_ntuple_{idx}.root") as f:
        f["Events"] = histo

print("Original events: ", original_events)
print("Cut events: ", cut_events)

os.system("hadd -f kit_ntuple_wmcr_MET_Jet.root kit_wmcr/kit_ntuple_*.root")
#os.system("hadd -f kit_ntuple_gcr_v4.root kit_gcr/kit_ntuple_*.root")
