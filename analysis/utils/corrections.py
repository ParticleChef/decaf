#! /usr/bin/env python
import correctionlib
from correctionlib import convert
import os
import awkward as ak

import numpy as np
from coffea import lookup_tools, jetmet_tools, util
from coffea.lookup_tools import extractor, dense_lookup
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory

import uproot
from coffea.util import save, load
import json

import hist

###
# MET trigger efficiency SFs, 2017/18 from monojet. Depends on recoil.
###

def get_met_trig_weight(year, met):
    met_trig_hists = {
        '2016postVFP': "data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root:hden_monojet_recoil_clone_passed",
        '2016preVFP': "data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root:hden_monojet_recoil_clone_passed",
        '2017': "data/trigger_eff/met_trigger_sf.root:120pfht_hltmu_1m_2017",
        '2018': "data/trigger_eff/met_trigger_sf.root:120pfht_hltmu_1m_2018"
    }
    corr = convert.from_uproot_THx(met_trig_hists[year])
    evaluator = corr.to_evaluator()

    met  = ak.fill_none(met, 0.)
    met  = ak.where((met>950.), ak.full_like(met,950.), met)

    weight = ak.where(
        ~np.isnan(ak.fill_none(met, np.nan)),
        evaluator.evaluate(met),
        ak.zeros_like(met)
    )
    return weight

####
# Electron ID scale factor
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_ele_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["UL-Electron-ID-SF"].evaluate(year, "sf", "Loose", flateta, flatpt)

    return ak.unflatten(weight, counts=counts)

def get_ele_tight_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["UL-Electron-ID-SF"].evaluate(year, "sf", "Tight", flateta, flatpt)
    
    return ak.unflatten(weight, counts=counts)


####
# Electron Trigger weight
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
# Copy from previous correctionsUL.py file
####

def get_ele_trig_weight(year, eta, pt):
    ele_trig_hists = {
        '2016postVFP': "data/ElectronTrigEff/egammaEffi.txt_EGM2D-2016postVFP.root:EGamma_SF2D",
        '2016preVFP' : "data/ElectronTrigEff/egammaEffi.txt_EGM2D-2016preVFP.root:EGamma_SF2D",
        '2017': "data/ElectronTrigEff/egammaEffi.txt_EGM2D-2017.root:EGamma_SF2D",#monojet measurement for the combined trigger path
        '2018': "data/ElectronTrigEff/egammaEffi.txt_EGM2D-2018.root:EGamma_SF2D" #approved by egamma group: https://indico.cern.ch/event/924522/
    }
    corr = convert.from_uproot_THx(ele_trig_hists[year])
    evaluator = corr.to_evaluator()

    eta = ak.fill_none(eta, 0.)
    pt = ak.fill_none(pt, 40.)
    pt  = ak.where((pt>250.),ak.full_like(pt,250.),pt)

    weight = ak.where(
        ~np.isnan(ak.fill_none(pt, np.nan)),
        evaluator.evaluate(eta, pt),
        ak.zeros_like(pt)
    )
    return weight

####
# Electron Reco scale factor
# root files: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
# Code Copy from previous correctionsUL.py file
####

def get_ele_reco_sf_below20(year, eta, pt):
    ele_reco_files_below20 = {
        '2016postVFP': "data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root:EGamma_SF2D",
        '2016preVFP': "data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root:EGamma_SF2D",
        '2017': "data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2017.root:EGamma_SF2D",
        '2018': "data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root:EGamma_SF2D"
    }

    corr = convert.from_uproot_THx(ele_reco_files_below20[year])
    evaluator = corr.to_evaluator()
    
    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    eta = ak.where((eta<2.399), ak.full_like(eta,-2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    pt = ak.where((pt>19.99), ak.full_like(pt,19.99), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator.evaluate(flateta, flatpt)
    return ak.unflatten(weight, counts=counts)
    #get_ele_reco_err_below20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.variances() ** 0.5, ele_reco_hist.axes)


def get_ele_reco_sf_above20(year, eta, pt):
    ele_reco_files_above20 = {
        '2016postVFP': "data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root:EGamma_SF2D",
        '2016preVFP': "data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root:EGamma_SF2D",
        '2017': "data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root:EGamma_SF2D",
        '2018': "data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root:EGamma_SF2D"
    }
    
    corr = convert.from_uproot_THx(ele_reco_files_above20[year])
    evaluator = corr.to_evaluator()
    
    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    eta = ak.where((eta<2.399), ak.full_like(eta,-2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<20.), ak.full_like(pt,20.), pt)
    pt = ak.where((pt>499.99), ak.full_like(pt,499.99), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator.evaluate(flateta, flatpt)
    return ak.unflatten(weight, counts=counts)
    #get_ele_reco_err_above20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.variances() ** 0.05, ele_reco_hist.axes)
    


####
# Photon ID scale factor
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_pho_tight_id_sf(year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/photon.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<20.),ak.full_like(pt,20.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["UL-Photon-ID-SF"].evaluate(year, "sf", "Tight", flateta, flatpt)

    return ak.unflatten(weight, counts=counts)


def get_pho_loose_id_sf(year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/photon.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt  = ak.where((pt<20.),ak.full_like(pt,20.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["UL-Photon-ID-SF"].evaluate(year, "sf", "Loose", flateta, flatpt)

    return ak.unflatten(weight, counts=counts)


#### 
# Photon CSEV sf 
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# 

#def get_pho_csev_sf(year, eta, pt):
#    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/photon.json.gz')
#
#    flateta, counts = ak.flatten(eta), ak.num(eta)
#    flatpt = ak.flatten(pt)
#    weight = evaluator["UL-Photon-CSEV-SF"].evaluate(year, "sf", "Tight", flateta, flatpt)
#
#    return ak.unflatten(weight, counts=counts)


####
# Photon Trigger weight
# Copy from previous decaf version
####

def get_pho_trig_weight(year, pt):
    pho_trig_files = {
        '2016postVFP': "data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root:hden_photonpt_clone_passed",
        '2016preVFP': "data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root:hden_photonpt_clone_passed",
        "2017": "data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root:hden_photonpt_clone_passed",
        "2018": "data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root:hden_photonpt_clone_passed"
    }

    corr = convert.from_uproot_THx(pho_trig_files[year])
    evaluator = corr.to_evaluator()

    return evaluator.evaluate(eta, pt)



# https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL

####
# Muon ID scale factor
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018n?topic=MuonUL2018
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/MUO
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_mu_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'_UL/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    if year == '2018':
        weight = evaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(year+'_UL', flateta, flatpt, "sf")
    else:
        weight = evaluator["NUM_LooseID_DEN_genTracks"].evaluate(year+'_UL', flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)

def get_mu_tight_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'_UL/muon_Z.json.gz')
    
    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    if year == '2018':
        weight = evaluator["NUM_TightID_DEN_TrackerMuons"].evaluate(year+'_UL', flateta, flatpt, "sf")
    else:
        weight = evaluator["NUM_TightID_DEN_genTracks"].evaluate(year+'_UL', flateta, flatpt, "sf")
    
    return ak.unflatten(weight, counts=counts)




####
# Muon Iso scale factor
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018n?topic=MuonUL2018
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/MUO
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_mu_loose_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'_UL/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_LooseRelIso_DEN_LooseID"].evaluate(year+'_UL', flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)

def get_mu_tight_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'_UL/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_TightRelIso_DEN_TightIDandIPCut"].evaluate(year+'_UL', flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)

###
# Muon scale and resolution (i.e. Rochester)
# https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
###

tag = 'roccor.Run2.v5'
get_mu_rochester_sf = {}
for year in ['2016postVFP', '2016preVFP', '2017','2018']:
    if '2016postVFP' in year: 
        fname = f'data/{tag}/RoccoR2016bUL.txt'
    elif '2016preVFP' in year:  
        fname = f'data/{tag}/RoccoR2016aUL.txt'
    else:
        fname = f'data/{tag}/RoccoR{year}UL.txt'
    sfs = lookup_tools.txt_converters.convert_rochester_file(fname,loaduncs=True)
    get_mu_rochester_sf[year] = lookup_tools.rochester_lookup.rochester_lookup(sfs)


####
# PU weight
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/LUM
####
#trueint = events.Pileup.nTrueInt
def get_pu_weight(year, trueint):
    correction = {'2018': 'Collisions18_UltraLegacy_goldenJSON',
                  '2017': 'Collisions17_UltraLegacy_goldenJSON',
                  '2016preVFP': 'Collisions16_UltraLegacy_goldenJSON',
                  '2016postVFP':'Collisions16_UltraLegacy_goldenJSON'}
    evaluator = correctionlib.CorrectionSet.from_file('data/PUweight/'+year+'_UL/puWeights.json.gz')
    weight = evaluator[correction[year]].evaluate(trueint, 'nominal')

    return weight


####
# XY MET Correction
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#xy_Shift_Correction_MET_phi_modu
####

# https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/
# correction_labels = ["metphicorr_pfmet_mc", "metphicorr_puppimet_mc", "metphicorr_pfmet_data", "metphicorr_puppimet_data"]

def XY_MET_Correction(year, npv, run, pt, phi, isData):
    
    npv = ak.where((npv>200),ak.full_like(npv,200),npv)
    pt  = ak.where((pt>1000.),ak.full_like(pt,1000.),pt)

    evaluator = correctionlib.CorrectionSet.from_file('data/JetMETCorr/'+year+'_UL/met.json.gz')

    if isData:
        corrected_pt = evaluator['pt_metphicorr_pfmet_data'].evaluate(pt,phi,npv,run)
        corrected_phi = evaluator['phi_metphicorr_pfmet_data'].evaluate(pt,phi,npv,run)

    if not isData:
        corrected_pt = evaluator['pt_metphicorr_pfmet_mc'].evaluate(pt,phi,npv,run)
        corrected_phi = evaluator['phi_metphicorr_pfmet_mc'].evaluate(pt,phi,npv,run)

    return corrected_pt, corrected_phi


####
# Jet
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
####


###
# V+jets NLO k-factors
# Only use nlo ewk sf
###

#nlo_ewk_hists = {
#    'dy': uproot.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_ewk"],
#    'w': uproot.open("data/vjets_SFs/merged_kfactors_wjets.root")["kfactor_monojet_ewk"],
#    'z': uproot.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_ewk"],
#    'a': uproot.open("data/vjets_SFs/merged_kfactors_gjets.root")["kfactor_monojet_ewk"]
#}    
#for p in ['dy','w','z','a']:
#    get_nlo_ewk_weight[p] = lookup_tools.dense_lookup.dense_lookup(nlo_ewk_hists[p].values(), nlo_ewk_hists[p].axes)
nlo_ewk_hists = {
        'dy': ["* * data/vjets_SFs/merged_kfactors_zjets.root"],
        'w': ["* * data/vjets_SFs/merged_kfactors_wjets.root"],
        'z': ["* * data/vjets_SFs/merged_kfactors_zjets.root"],
        'a': ["* * data/vjets_SFs/merged_kfactors_gjets.root"],
}    
get_nlo_ewk_weight = {}
for p in ['dy','w','z','a']:
    print(nlo_ewk_hists[p])
    ext = extractor()
    ext.add_weight_sets(nlo_ewk_hists[p])
    ext.finalize()
    get_nlo_ewk_weight[p] = ext.make_evaluator()["kfactor_monojet_ewk"]

###
# V+jets NNLO weights
# The schema is process_NNLO_NLO_QCD1QCD2QCD3_EW1EW2EW3_MIX, where 'n' stands for 'nominal', 'u' for 'up', and 'd' for 'down'
###

histname={
    'dy': 'eej_NNLO_NLO_',
    'w':  'evj_NNLO_NLO_',
    'z': 'vvj_NNLO_NLO_',
    'a': 'aj_NNLO_NLO_'
}
correlated_variations = {
    'cen':    'nnn_nnn_n',
    'qcd1up': 'unn_nnn_n',
    'qcd1do': 'dnn_nnn_n',
    'qcd2up': 'nun_nnn_n',
    'qcd2do': 'ndn_nnn_n',
    'qcd3up': 'nnu_nnn_n',
    'qcd3do': 'nnd_nnn_n',
    'ew1up' : 'nnn_unn_n',
    'ew1do' : 'nnn_dnn_n',
    'mixup' : 'nnn_nnn_u',
    'mixdo' : 'nnn_nnn_d',
    'muFup' : 'nnn_nnn_n_Weight_scale_variation_muR_1p0_muF_2p0',
    'muFdo' : 'nnn_nnn_n_Weight_scale_variation_muR_1p0_muF_0p5',
    'muRup' : 'nnn_nnn_n_Weight_scale_variation_muR_2p0_muF_1p0',
    'muRdo' : 'nnn_nnn_n_Weight_scale_variation_muR_0p5_muF_1p0'
}
uncorrelated_variations = {
    'dy': {
        'ew2Gup': 'nnn_nnn_n',
        'ew2Gdo': 'nnn_nnn_n',
        'ew2Wup': 'nnn_nnn_n',
        'ew2Wdo': 'nnn_nnn_n',
        'ew2Zup': 'nnn_nun_n',
        'ew2Zdo': 'nnn_ndn_n',
        'ew3Gup': 'nnn_nnn_n',
        'ew3Gdo': 'nnn_nnn_n',
        'ew3Wup': 'nnn_nnn_n',
        'ew3Wdo': 'nnn_nnn_n',
        'ew3Zup': 'nnn_nnu_n',
        'ew3Zdo': 'nnn_nnd_n'
    },
    'w': {
        'ew2Gup': 'nnn_nnn_n',
        'ew2Gdo': 'nnn_nnn_n',
        'ew2Wup': 'nnn_nun_n',
        'ew2Wdo': 'nnn_ndn_n',
        'ew2Zup': 'nnn_nnn_n',
        'ew2Zdo': 'nnn_nnn_n',
        'ew3Gup': 'nnn_nnn_n',
        'ew3Gdo': 'nnn_nnn_n',
        'ew3Wup': 'nnn_nnu_n',
        'ew3Wdo': 'nnn_nnd_n',
        'ew3Zup': 'nnn_nnn_n',
        'ew3Zdo': 'nnn_nnn_n'
    },
    'z': {
        'ew2Gup': 'nnn_nnn_n',
        'ew2Gdo': 'nnn_nnn_n',
        'ew2Wup': 'nnn_nnn_n',
        'ew2Wdo': 'nnn_nnn_n',
        'ew2Zup': 'nnn_nun_n',
        'ew2Zdo': 'nnn_ndn_n',
        'ew3Gup': 'nnn_nnn_n',
        'ew3Gdo': 'nnn_nnn_n',
        'ew3Wup': 'nnn_nnn_n',
        'ew3Wdo': 'nnn_nnn_n',
        'ew3Zup': 'nnn_nnu_n',
        'ew3Zdo': 'nnn_nnd_n'
    },
    'a': {
        'ew2Gup': 'nnn_nun_n',
        'ew2Gdo': 'nnn_ndn_n',
        'ew2Wup': 'nnn_nnn_n',
        'ew2Wdo': 'nnn_nnn_n',
        'ew2Zup': 'nnn_nnn_n',
        'ew2Zdo': 'nnn_nnn_n',
        'ew3Gup': 'nnn_nnu_n',
        'ew3Gdo': 'nnn_nnd_n',
        'ew3Wup': 'nnn_nnn_n',
        'ew3Wdo': 'nnn_nnn_n',
        'ew3Zup': 'nnn_nnn_n',
        'ew3Zdo': 'nnn_nnn_n'
    }
}
get_nnlo_nlo_weight = {}
for year in ['2016postVFP', '2016preVFP', '2017','2018']:
    get_nnlo_nlo_weight[year] = {}
    if '2016' in year:
        nnlo_file = {
            'dy': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_eej_madgraph_2016.root"),
            'w': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_evj_madgraph_2016.root"),
            'z': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_vvj_madgraph_2016.root"),
            'a': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_aj_madgraph_2016.root")
        }
    else:
        nlo_file = {
            'dy': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_eej_madgraph_"+year+".root"),
            'w': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_evj_madgraph_"+year+".root"),
            'z': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_vvj_madgraph_"+year+".root"),
            'a': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_aj_madgraph_"+year+".root")
        }
    for p in ['dy','w','z','a']:
        get_nnlo_nlo_weight[year][p] = {}
        for cv in correlated_variations:
            histo=nnlo_file[p][histname[p]+correlated_variations[cv]]
            get_nnlo_nlo_weight[year][p][cv]=lookup_tools.dense_lookup.dense_lookup(histo.values(), histo.axes)
        for uv in uncorrelated_variations[p]:
            histo=nnlo_file[p][histname[p]+uncorrelated_variations[p][uv]]
            get_nnlo_nlo_weight[year][p][uv]=lookup_tools.dense_lookup.dense_lookup(histo.values(), histo.axes)



def get_ttbar_weight(pt):
    return np.exp(0.0615 - 0.0005 * np.clip(pt, 0, 800))


# Soft drop mass correction updated for UL. Copied from:
# https://github.com/jennetd/hbb-coffea/blob/master/boostedhiggs/corrections.py
# Renamed copy of corrected_msoftdrop() function

msdcorr = correctionlib.CorrectionSet.from_file('data/msdcorr.json')

def get_msd_corr(fatjets):
    msdraw = np.sqrt(
        np.maximum(
            0.0,
            (fatjets.subjets * (1 - fatjets.subjets.rawFactor)).sum().mass2,
        )
    )
    msdfjcorr = msdraw / (1 - fatjets.rawFactor)

    corr = msdcorr["msdfjcorr"].evaluate(
        np.array(ak.flatten(msdfjcorr / fatjets.pt)),
        np.array(ak.flatten(np.log(fatjets.pt))),
        np.array(ak.flatten(fatjets.eta)),
    )
    corr = ak.unflatten(corr, ak.num(fatjets))
    corrected_mass = msdfjcorr * corr

    return corrected_mass

from coffea.lookup_tools.correctionlib_wrapper import correctionlib_wrapper
from coffea.lookup_tools.dense_lookup import dense_lookup

class BTagCorrector:

    def __init__(self, tagger, year, workingpoint):
        self._year = year

        wp = {}
        wp['loose'] = 'L'
        wp['medium'] = 'M'
        wp['tight'] = 'T'
        self._wp = wp[workingpoint]

        btvjson = {}
        btvjson['deepflav'] = {
            'incl': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'_UL/btagging.json.gz')["deepJet_incl"],
            'comb': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'_UL/btagging.json.gz')["deepJet_comb"],
        }
        btvjson['deepcsv'] = {
            'incl': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'_UL/btagging.json.gz')["deepCSV_incl"],
            'comb': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'_UL/btagging.json.gz')["deepCSV_comb"],
        }
        self.sf = btvjson[tagger]

        files = {
            '2016preVFP': 'btageff2016preVFP.merged',
            '2016postVFP': 'btageff2016postVFP.merged',
            '2017': 'btageff2017.merged',
            '2018': 'btageff2018.merged',
        }
        filename = 'hists/'+files[year]
        btag_file = load(filename)
        for k in btag_file[tagger]:
            try:
                btag += btag_file[tagger][k]
            except:
                btag = btag_file[tagger][k]
        bpass = btag[{"wp": workingpoint, "btag": "pass"}].view()
        ball = btag[{"wp": workingpoint, "btag": sum}].view()
        ball[ball<=0.]=1.
        ratio = bpass / np.maximum(ball, 1.)
        nom = hist.Hist(*btag.axes[2:], data=ratio)
        nom.name = "ratios"  
        nom.label = "out"
        self.eff = convert.from_histogram(nom).to_evaluator()

    def btag_weight(self, pt, eta, flavor, istag):

        abseta = abs(eta)
        flateta, counts = ak.fill_none(ak.flatten(abseta), 0.), ak.num(abseta)

        pt = ak.where((pt>999.99), ak.full_like(pt,999.99), pt)
        flatpt =  ak.fill_none(ak.flatten(pt), 30.)

        flatflavor = ak.fill_none(ak.flatten(flavor), 0)
        
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1b_Event_reweighting_using_scale
        def P(eff):
            weight = ak.where(istag, eff, 1-eff)
            return ak.prod(weight, axis=1)

        '''
        Correction deepJet_comb has 5 inputs
        Input systematic (string): 
        Input working_point (string): L/M/T
        Input flavor (int): hadron flavor definition: 5=b, 4=c, 0=udsg
        Input abseta (real):
        Input pt (real):
        '''

        eff = ak.where(
            ~np.isnan(ak.fill_none(pt, np.nan)),
            ak.unflatten(self.eff.evaluate(flatflavor, flatpt, flateta), counts=counts),
            ak.zeros_like(pt)
        )
        sf_nom = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_up_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_down_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_up_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_down_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('down_uncorrelated',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('down_uncorrelated',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)        
            )
        )
        sf_light_up_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_down_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_up_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_down_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['incl'].evaluate('down_uncorrelated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        
        eff_data_nom  = ak.where(
            (sf_nom*eff>1.), 
            ak.ones_like(eff), 
            sf_nom*eff
        )
        eff_data_bc_up_correlated   = ak.where(
            (sf_bc_up_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_up_correlated*eff
        )
        eff_data_bc_down_correlated = ak.where(
            (sf_bc_down_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_down_correlated*eff
        )
        eff_data_bc_up_uncorrelated = ak.where(
            (sf_bc_up_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_up_uncorrelated*eff
        )
        eff_data_bc_down_uncorrelated = ak.where(
            (sf_bc_down_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_down_uncorrelated*eff
        )
        eff_data_light_up_correlated   = ak.where(
            (sf_light_up_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_up_correlated*eff
        )
        eff_data_light_down_correlated = ak.where(
            (sf_light_down_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_down_correlated*eff
        )
        eff_data_light_up_uncorrelated = ak.where(
            (sf_light_up_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_up_uncorrelated*eff
        )
        eff_data_light_down_uncorrelated = ak.where(
            (sf_light_down_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_down_uncorrelated*eff
        )
       
        nom = P(eff_data_nom)/P(eff)
        bc_up_correlated = P(eff_data_bc_up_correlated)/P(eff)
        bc_down_correlated = P(eff_data_bc_down_correlated)/P(eff)
        bc_up_uncorrelated = P(eff_data_bc_up_uncorrelated)/P(eff)
        bc_down_uncorrelated = P(eff_data_bc_down_uncorrelated)/P(eff)
        light_up_correlated = P(eff_data_light_up_correlated)/P(eff)
        light_down_correlated = P(eff_data_light_down_correlated)/P(eff)
        light_up_uncorrelated = P(eff_data_light_up_uncorrelated)/P(eff)
        light_down_uncorrelated = P(eff_data_light_down_uncorrelated)/P(eff)
        
        return np.nan_to_num(nom, nan=1.), \
        np.nan_to_num(bc_up_correlated, nan=1.), \
        np.nan_to_num(bc_down_correlated, nan=1.), \
        np.nan_to_num(bc_up_uncorrelated, nan=1.), \
        np.nan_to_num(bc_down_uncorrelated, nan=1.), \
        np.nan_to_num(light_up_correlated, nan=1.), \
        np.nan_to_num(light_down_correlated, nan=1.), \
        np.nan_to_num(light_up_uncorrelated, nan=1.), \
        np.nan_to_num(light_down_uncorrelated, nan=1.)

jec_name_map = {
    'JetPt': 'pt',
    'JetMass': 'mass',
    'JetEta': 'eta',
    'JetA': 'area',
    'ptGenJet': 'pt_gen',
    'ptRaw': 'pt_raw',
    'massRaw': 'mass_raw',
    'Rho': 'event_rho',
    'METpt': 'pt',
    'METphi': 'phi',
    'JetPhi': 'phi',
    'UnClusteredEnergyDeltaX': 'MetUnclustEnUpDeltaX',
    'UnClusteredEnergyDeltaY': 'MetUnclustEnUpDeltaY',
}

def jet_factory_factory(files):
    ext = extractor()
    directory='data/jerc'
    for filename in files:
        ext.add_weight_sets([f"* * {directory+'/'+filename}"])
    ext.finalize()
    jec_stack = JECStack(ext.make_evaluator())
    return CorrectedJetsFactory(jec_name_map, jec_stack)

jet_factory = {
    "2016preVFPmc": jet_factory_factory(
        files=[
            "Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
            "Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.junc.txt",
            "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
            "Summer20UL16APV_JRV3_MC_SF_AK4PFchs.jersf.txt",
        ]
    ),
    "2016preVFPmcNOJER": jet_factory_factory(
        files=[
            "Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.junc.txt",
        ]
    ),
    "2016postVFPmc": jet_factory_factory(
        files=[
            "Summer19UL16_V7_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL16_V7_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
            "Summer19UL16_V7_MC_Uncertainty_AK4PFchs.junc.txt",
            "Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
            "Summer20UL16_JRV3_MC_SF_AK4PFchs.jersf.txt",
        ]
    ),
    "2016postVFPmcNOJER": jet_factory_factory(
        files=[
            "Summer19UL16_V7_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL16_V7_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL16_V7_MC_Uncertainty_AK4PFchs.junc.txt",
        ]
    ),
    "2017mc": jet_factory_factory(
        files=[
            "Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
            "Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt",
            "Summer19UL17_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
            "Summer19UL17_JRV3_MC_SF_AK4PFchs.jersf.txt",
        ]
    ),
    "2017mcNOJER": jet_factory_factory(
        files=[
            "Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt",
        ]
    ),
    "2018mc": jet_factory_factory(
        files=[
            "Summer19UL18_V5_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL18_V5_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
            "Summer19UL18_V5_MC_Uncertainty_AK4PFchs.junc.txt",
            "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
            "Summer19UL18_JRV2_MC_SF_AK4PFchs.jersf.txt",
        ]
    ),
    "2018mcNOJER": jet_factory_factory(
        files=[
            "Summer19UL18_V5_MC_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL18_V5_MC_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL18_V5_MC_Uncertainty_AK4PFchs.junc.txt",
        ]
    ),
}

subjet_factory = {
    "2016preVFPmc": jet_factory_factory(
        files=[
            "Summer19UL16APV_V7_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_UncertaintySources_AK4PFPuppi.junc.txt",
            "Summer19UL16APV_V7_MC_Uncertainty_AK4PFPuppi.junc.txt",
            "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFPuppi.jr.txt",
            "Summer20UL16APV_JRV3_MC_SF_AK4PFPuppi.jersf.txt",
        ]
    ),
    "2016preVFPmcNOJER": jet_factory_factory(
        files=[
            "Summer19UL16APV_V7_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_Uncertainty_AK4PFPuppi.junc.txt",
        ]
    ),
    "2016postVFPmc": jet_factory_factory(
        files=[
            "Summer19UL16_V7_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_UncertaintySources_AK4PFPuppi.junc.txt",
            "Summer19UL16_V7_MC_Uncertainty_AK4PFPuppi.junc.txt",
            "Summer20UL16_JRV3_MC_PtResolution_AK4PFPuppi.jr.txt",
            "Summer20UL16_JRV3_MC_SF_AK4PFPuppi.jersf.txt",
        ]
    ),
    "2016postVFPmcNOJER": jet_factory_factory(
        files=[
            "Summer19UL16_V7_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_Uncertainty_AK4PFPuppi.junc.txt",
        ]
    ),
    "2017mc": jet_factory_factory(
        files=[
            "Summer19UL17_V5_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_UncertaintySources_AK4PFPuppi.junc.txt",
            "Summer19UL17_V5_MC_Uncertainty_AK4PFPuppi.junc.txt",
            "Summer19UL17_JRV3_MC_PtResolution_AK4PFPuppi.jr.txt",
            "Summer19UL17_JRV3_MC_SF_AK4PFPuppi.jersf.txt",
        ]
    ),
    "2017mcNOJER": jet_factory_factory(
        files=[
            "Summer19UL17_V5_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_Uncertainty_AK4PFPuppi.junc.txt",
        ]
    ),
    "2018mc": jet_factory_factory(
        files=[
            "Summer19UL18_V5_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_UncertaintySources_AK4PFPuppi.junc.txt",
            "Summer19UL18_V5_MC_Uncertainty_AK4PFPuppi.junc.txt",
            "Summer19UL18_JRV2_MC_PtResolution_AK4PFPuppi.jr.txt",
            "Summer19UL18_JRV2_MC_SF_AK4PFPuppi.jersf.txt",
        ]
    ),
    "2018mcNOJER": jet_factory_factory(
        files=[
            "Summer19UL18_V5_MC_L1FastJet_AK4PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_L2Relative_AK4PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_Uncertainty_AK4PFPuppi.junc.txt",
        ]
    ),
}

fatjet_factory = {
    "2016preVFPmc": jet_factory_factory(
        files=[
            "Summer19UL16APV_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_UncertaintySources_AK8PFPuppi.junc.txt",
            "Summer19UL16APV_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
            "Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
            "Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
        ]
    ),
    "2016preVFPmcNOJER": jet_factory_factory(
        files=[
            "Summer19UL16APV_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL16APV_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
        ]
    ),
    "2016postVFPmc": jet_factory_factory(
        files=[
            "Summer19UL16_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_UncertaintySources_AK8PFPuppi.junc.txt",
            "Summer19UL16_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
            "Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
            "Summer20UL16_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
        ]
    ),
    "2016postVFPmcNOJER": jet_factory_factory(
        files=[
            "Summer19UL16_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL16_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
        ]
    ),
    "2017mc": jet_factory_factory(
        files=[
            "Summer19UL17_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_UncertaintySources_AK8PFPuppi.junc.txt",
            "Summer19UL17_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
            "Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
            "Summer19UL17_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
        ]
    ),
    "2017mcNOJER": jet_factory_factory(
        files=[
            "Summer19UL17_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL17_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
        ]
    ),
    "2018mc": jet_factory_factory(
        files=[
            "Summer19UL18_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_UncertaintySources_AK8PFPuppi.junc.txt",
            "Summer19UL18_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
            "Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.jr.txt",
            "Summer19UL18_JRV2_MC_SF_AK8PFPuppi.jersf.txt",
        ]
    ),
    "2018mcNOJER": jet_factory_factory(
        files=[
            "Summer19UL18_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
            "Summer19UL18_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
        ]
    ),
}
met_factory = CorrectedMETFactory(jec_name_map)


corrections = {}
corrections = {
    'get_met_trig_weight':      get_met_trig_weight,
    'get_ele_loose_id_sf':      get_ele_loose_id_sf,
    'get_ele_tight_id_sf':      get_ele_tight_id_sf,
    'get_ele_trig_weight':      get_ele_trig_weight,
    'get_ele_reco_sf_below20':  get_ele_reco_sf_below20,
    #'get_ele_reco_err_below20': get_ele_reco_err_below20,
    'get_ele_reco_sf_above20':  get_ele_reco_sf_above20,
    #'get_ele_reco_err_above20': get_ele_reco_err_above20,
    'get_pho_loose_id_sf':      get_pho_loose_id_sf,
    'get_pho_tight_id_sf':      get_pho_tight_id_sf,
    'get_pho_trig_weight':      get_pho_trig_weight,
    'get_mu_loose_id_sf':       get_mu_loose_id_sf,
    'get_mu_tight_id_sf':       get_mu_tight_id_sf,
    'get_mu_loose_iso_sf':      get_mu_loose_iso_sf,
    'get_mu_tight_iso_sf':      get_mu_tight_iso_sf,
    'get_met_xy_correction':    XY_MET_Correction,
    'get_pu_weight':            get_pu_weight,
    'get_nlo_ewk_weight':       get_nlo_ewk_weight,
    'get_nnlo_nlo_weight':      get_nnlo_nlo_weight,
    'get_ttbar_weight':         get_ttbar_weight,
    'get_msd_corr':             get_msd_corr,
    'get_btag_weight':          BTagCorrector,
    'get_mu_rochester_sf':      get_mu_rochester_sf,
    'jet_factory':              jet_factory,
    'subjet_factory':           subjet_factory,
    'fatjet_factory':           fatjet_factory,
    'met_factory':              met_factory
}


save(corrections, 'data/corrections.coffea')
