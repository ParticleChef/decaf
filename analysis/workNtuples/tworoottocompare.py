import uproot
import numpy as np
import awkward as ak
import json
import gzip
import tqdm
import vector
from multiprocessing import Pool, Manager
from ntuplelist import *
from coffea.util import load, save
from present_from_taiwoo import pid
from custom_dr import *
corrections = load('/home/jhong/monotop/CMSSW_11_3_4/src/decaf/analysis/data/corrections.coffea')

year = '2018'
isRealsample = True
sample_name = 'A' ## A, B, C, D, ALL
branches_to_load = ['N_HEMJets_nom', 'RawJet_Eta_nom', 'RawJet_Phi_nom',
                    'AK15Jet_Pt_nom', 'AK15Jet_chHEF_nom', 'AK15Jet_neHEF_nom', 'N_AK15Jets_nom', 'N_AK15SoftDropJets_nom',
                    'N_LooseElectrons', 'N_LooseMuons', 'N_LoosePhotons', 'N_TightPhotons',
                    'Jets_outside_lead_AK15Jet_taggedL_nom', 'AK15Jet_particleNet_TvsQCD_nom',
                    'TightMuon_Pt', 'TightMuon_Eta', 'TightMuon_Phi', 'LooseMuon_Pt', 'LooseMuon_Eta', 'LooseMuon_Phi'
                    ]

def GetEventID(f,Lumi, Run, Event):
    lumi = f["Events/"+Lumi].array()
    run = f["Events/"+Run].array()
    event = f["Events/"+Event].array()
    return lumi, run, event

Moritz = uproot.open("/home/jhong/monotop/CMSSW_11_3_4/src/decaf/analysis/workNtuples/kit_ntuple_wmcr_MET.root")
lumi, run, event = GetEventID(Moritz, "Lumi", "Run", "Event")

exr = event / run
print(len(exr), len(np.unique(exr))) 

metadata = '/home/jhong/monotop/CMSSW_11_3_4/src/decaf/analysis/metadata/KNUv1_UL_ALL2018_v4.json.gz'
flist = []
with gzip.open(metadata, 'r') as f:
    metadata = json.load(f)
    target = 'MET'#+sample_name
    for key in metadata.keys():
        if target in key:
            print(key)
            samplename = key.split('_')[0]
            files = metadata[key]['files']
        else:
            continue
        for file in files:
            flist.append(file)

#with uproot.open(sample[sample_name]) as myfile:

def process_file(args):
    file = args
    nEvents = 0
    histo = {}
    get_met_xy_correction = corrections['get_met_xy_correction']
    with uproot.open(file) as f:
        outname = file.split('/')[-1].split('.')[0].split('_')[-1]
        samplename = file.split('/')[-2]
        tr = f['Events']
        run = tr['run'].array()
        event = tr['event'].array()
        npv = tr['PV_npvsGood'].array()
        electron = tr.arrays(["Electron_pt", "Electron_eta", "Electron_phi", "Electron_cutBased","Electron_mass","Electron_charge","Electron_deltaEtaSC"], how="zip")['Electron']
        photon = tr.arrays(["Photon_pt", "Photon_eta", "Photon_phi", "Photon_cutBased", "Photon_mass", "Photon_charge","Photon_electronVeto","Photon_isScEtaEB"], how="zip")['Photon']
        muon = tr.arrays(["Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge", "Muon_looseId", "Muon_tightId", "Muon_pfRelIso04_all", "Muon_pfIsoId", "Muon_isTracker", "Muon_isPFcand","Muon_isGlobal"], how="zip")['Muon']
        ak15 = tr.arrays(["AK15PFPuppi_Jet_pt", "AK15PFPuppi_Jet_eta", "AK15PFPuppi_Jet_phi", "AK15PFPuppi_Jet_mass", "AK15PFPuppi_Jet_jetId", "AK15PFPuppi_Jet_neHEF", "AK15PFPuppi_Jet_chHEF",
                          "AK15PFPuppi_Jet_particleNetAK15_QCDbb","AK15PFPuppi_Jet_particleNetAK15_QCDcc","AK15PFPuppi_Jet_particleNetAK15_QCDb","AK15PFPuppi_Jet_particleNetAK15_QCDc","AK15PFPuppi_Jet_particleNetAK15_QCDothers","AK15PFPuppi_Jet_particleNetAK15_Tbqq","AK15PFPuppi_Jet_particleNetAK15_Tbcq"], how="zip")['AK15PFPuppi_Jet']
        ak4 = tr.arrays(["Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId", "Jet_btagDeepFlavB","Jet_btagDeepB","Jet_puId", "Jet_neHEF", "Jet_chHEF"], how="zip")['Jet']
        met = tr.arrays(["MET_pt", "MET_phi"], how="zip")
        METXY_pt, METXY_phi = get_met_xy_correction('2018', npv, run, met['MET_pt'], met['MET_phi'], isData=True)
        flat_variables = tr.arrays(["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight","HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_IsoMu24","Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter","Flag_BadPFMuonDzFilter","Flag_eeBadScFilter","Flag_ecalBadCalibFilter","PV_npvsGood","run", "luminosityBlock","event"], how="zip")
        flat_variables['METXY_pt'] = METXY_pt
        flat_variables['METXY_phi'] = METXY_phi
        flat_variables['MET_pt'] = met["MET_pt"]
        flat_variables['MET_phi'] = met["MET_phi"]
        single_photon = tr['HLT_Photon200'].array()
        single_electron = tr['HLT_Ele32_WPTight_Gsf'].array()
        PFMETNoMu120_PFHT60 = tr['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60'].array()
        PFMETNoMu120 = tr['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'].array()
        # run / event values should be matched with "exr" from Moritz
        myexr = event / run
        e_loose = pid.isLooseElectron(electron)
        e_tight = pid.isTightElectron(electron)
        a_loose = pid.isLoosePhoton(photon)
        a_tight = pid.isTightPhoton(photon)
        m_loose = pid.isLooseMuon(muon)
        m_tight = pid.isTightMuon(muon)
        j_good = pid.isGoodAK4(ak4)
        j_hem = pid.isHEMJet(ak4)
        fj_good = pid.isGoodAK15(ak15)
        n_j = ak.num(j_good)
        n_fj = ak.num(fj_good)
        n_eloose = ak.num(e_loose)
        n_etight = ak.num(e_tight)
        n_aloose = ak.num(a_loose)
        n_atight = ak.num(a_tight)
        n_mloose = ak.num(m_loose)
        n_mtight = ak.num(m_tight)
        n_HEM = ak.num(j_hem)

        filter_2018 = tr['Flag_goodVertices'].array() & tr['Flag_globalSuperTightHalo2016Filter'].array() & tr['Flag_HBHENoiseFilter'].array() & tr['Flag_HBHENoiseIsoFilter'].array() & tr['Flag_EcalDeadCellTriggerPrimitiveFilter'].array() & tr['Flag_BadPFMuonFilter'].array() & tr['Flag_BadPFMuonDzFilter'].array() & tr['Flag_eeBadScFilter'].array() & tr['Flag_ecalBadCalibFilter'].array()
        single_electron_trigger = (single_photon | single_electron)
        met_trigger = (PFMETNoMu120_PFHT60 | PFMETNoMu120)

        cr_mask = (n_eloose == 0) & (n_aloose == 0) & (n_atight == 0) & (n_mtight == 1) & (n_mloose == 1)
        mask = np.isin(myexr, exr) & filter_2018 & met_trigger & (n_HEM == 0) & (n_fj > 0) #& cr_mask 

        flat_variables = flat_variables[mask]
        e_loose, e_tight, a_loose, a_tight = e_loose[mask], e_tight[mask], a_loose[mask], a_tight[mask]
        j_good,  fj_good = j_good[mask],  fj_good[mask]
        m_loose, m_tight = m_loose[mask], m_tight[mask]
        a_raw, e_raw = photon[mask], electron[mask]
        m_raw = muon[mask]

        recoil_mask = (fj_good.pt[:,0] > 250)
        if len(flat_variables['MET_pt']) == 0:
            return 0

        m_vector = ak.zip({
            "pt": m_loose.pt[:,0],
            "eta": m_loose.eta[:,0],
            "phi": m_loose.phi[:,0],
            "mass": m_loose.mass[:,0],
            },
            with_name="LorentzVector",
            behavior=ak.behavior,
        )
        e_vector = ak.zip({
            "pt": e_loose.pt,
            "eta": e_loose.eta,
            "phi": e_loose.phi,
            "mass": e_loose.mass,
            },
            with_name="LorentzVector",
            behavior=ak.behavior,
        )
        fj_vector = ak.zip({
            "pt": fj_good.pt[:,0],
            "eta": fj_good.eta[:,0],
            "phi": fj_good.phi[:,0],
            "mass": fj_good.mass[:,0],
            },
            with_name="LorentzVector",
            behavior=ak.behavior,
        )
        dr_fjm = delta_r(fj_vector, m_vector)
        #print("dr: ", ak.type(dr_fjm))
        #print("fj: ", ak.type(fj_good))
        fj_good['isclean'] = (dr_fjm > 1.5)
        #print(fj_good.isclean)

        fj_clean = fj_good[fj_good.isclean]

        histo[''] = flat_variables
        histo['RawElectron'] = e_raw
        histo['RawPhoton'] = a_raw
        histo['RawMuon'] = m_raw
        histo['LooseElectron'] = e_loose
        histo['TightElectron'] = e_tight
        histo['LoosePhoton'] = a_loose
        histo['TightPhoton'] = a_tight
        histo['LooseMuon'] = m_loose
        histo['TightMuon'] = m_tight
        histo['Muon'] = m_loose
        histo['Jet'] = j_good
        histo['FJetGood'] = fj_good
        histo['FJet'] = fj_clean
        histo['dR_fj_m'] = dr_fjm
#        histo['AK15PFPuppi_Jet'] = fj_clean
        with uproot.recreate("/home/jhong/monotop/CMSSW_11_3_4/src/decaf/analysis/workNtuples/knu_wmcr/"+samplename+"_"+outname+".root") as f:
            f["Events"] = histo
    return len(flat_variables['MET_phi'])        

def main(flist):
    with Pool() as pool:
        results = list(tqdm.tqdm(pool.imap(
            process_file, [(file) for file in flist]), total=len(flist)))
        # sum all the results
        total = sum(results)
        print("Total events:", total)

if __name__ == "__main__":
    main(flist)

'''
for mytree in uproot.iterate(
    {file: "Events" for file in sample[sample_name]},
    branches=branches_to_load,
    library="ak",
    step_size= 120000 #"200 MB"
):
#    mytree = myfile['Events']
    Evt_ID = mytree['Evt_ID']
    Evt_Run = mytree['Evt_Run']

    TightMuon_Pt  = mytree['TightMuon_Pt']
    TightMuon_Eta = mytree['TightMuon_Eta']
    TightMuon_Phi = mytree['TightMuon_Phi']
    LooseMuon_Pt  = mytree['LooseMuon_Pt']
    LooseMuon_Eta = mytree['LooseMuon_Eta']
    LooseMuon_Phi = mytree['LooseMuon_Phi']

    onearray = ak.Array(np.ones_like(HLT_Photon200))
    noHEM = (RawJet_Eta_nom < -1.3) & (RawJet_Eta_nom > -3.0) & (RawJet_Phi_nom < -0.87) & (RawJet_Phi_nom > -1.57)
    n_total = n_total + ak.sum(onearray)
    
    cut1_HEM = ( N_HEMJets_nom == 0) & (ak.sum(noHEM,axis=1)==0)# & (
    n_HEM = n_HEM + ak.sum(onearray[cut1_HEM])

    cut2_ak15 = cut1_HEM & (ak.firsts(AK15Jet_Pt_nom) >= 250.) & (N_AK15Jets_nom >= 1) & (N_Jets_nom >= 1) & (ak.firsts(AK15Jet_chHEF_nom) > 0.1) & (ak.firsts(AK15Jet_neHEF_nom) < 0.8) & (N_AK15SoftDropJets_nom == N_AK15Jets_nom)
    n_ak15 = n_ak15 + ak.sum(onearray[cut2_ak15])
    
    cut3_recoil = cut2_ak15 & (Hadr_Recoil_MET_T1XY_Pt_nom >= 350.) & (ak.firsts(DeltaPhi_AK15Jet_Hadr_Recoil_MET_T1XY_nom) > 1.5) & ak.all(DeltaPhi_AK4Jet_Hadr_Recoil_MET_T1XY_nom > 0.5, axis=1)
    n_recoil = n_recoil + ak.sum(onearray[cut3_recoil])
    
    cut4_trig = cut3_recoil & (HLT_Photon200 == 1)
    n_trigger = n_trigger + ak.sum(onearray[cut4_trig]) 
    
    cut_gcr = cut4_trig & (N_LooseElectrons == 0) & (N_LooseMuons == 0) & (N_LoosePhotons == 1) & (N_TightPhotons == 1) & (ak.sum(Jets_outside_lead_AK15Jet_taggedL_nom, axis=1) == 0)
    n_gcr = n_gcr + ak.sum(onearray[cut_gcr])

    cut_pass = (AK15Jet_particleNet_TvsQCD_nom >= 0.26)
    cut_fail = (AK15Jet_particleNet_TvsQCD_nom < 0.26)

def print_cut():
    print('-----------------------------')
    print(n_total)
    print(n_HEM)
    print(n_ak15)
    print(n_recoil)
    print(n_trigger)
    print(n_gcr)
    #print(cut_pass)
    #print(cut_fail)


print_cut()
'''

