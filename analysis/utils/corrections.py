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

####
# PU weight
# https://cms-analysis-corrections.docs.cern.ch/corrections_era/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/LUM/latest/#puweights_cjsongz
####
#trueint = events.Pileup.nTrueInt
def get_pu_weight(year, trueint):
    correction = {
        '2022pre' : 'Collisions2022_355100_357900_eraBCD_GoldenJson',
        '2022post': 'Collisions2022_359022_362760_eraEFG_GoldenJson',
        '2023pre' : 'Collisions2023_366403_369802_eraBC_GoldenJson',
        '2023post': 'Collisions2023_369803_370790_eraD_GoldenJson',
        #'2024': 'Collisions2024_378981_386951_GoldenJson',
        '2024': 'Collisions24_BCDEFGHI_goldenJSON',
    }
    evaluator = correctionlib.CorrectionSet.from_file('data/PUweight/'+year+'/puWeights.json.gz')
    weight = evaluator[correction[year]].evaluate(trueint, 'nominal')
    systup = evaluator[correction[year]].evaluate(trueint, 'up')
    systdown = evaluator[correction[year]].evaluate(trueint, 'down')

    return weight, systup, systdown

# MET XY Correction
def get_met_xy_correction(year, met_type, isData, met_pt, met_phi, npvGood):
    evaluator = correctionlib.CorrectionSet.from_file('data/JMESF/'+year+'/met_xy_correction.json.gz')
    if year == '2022pre':
        epoch = '2022'
    elif year == '2022post':
        epoch = '2022EE'
    else:
        pass
    if isData:
        dtmc = 'DATA'
    else:
        dtmc = 'MC'
    corr_pt = evaluator['met_xy_corrections'].evaluate('pt', met_type, epoch, dtmc, 'nom', met_pt, met_phi, npvGood)
    corr_phi = evaluator['met_xy_corrections'].evaluate('phi', met_type, epoch, dtmc, 'nom', met_pt, met_phi, npvGood)
    return corr_pt, corr_phi

def get_jec_correction(year, pt, eta, phi, rho, area, run, isData):
    evaluator = correctionlib.CorrectionSet.from_file('data/JMESF/'+year+'/jet_jerc.json.gz')
    counts = ak.num(pt)
    run, _ = ak.broadcast_arrays(run, pt)
    rho, _ = ak.broadcast_arrays(rho, pt)
    pt, eta, phi, rho, area, run = ak.flatten(pt), ak.flatten(eta), ak.flatten(phi), ak.flatten(rho), ak.flatten(area), ak.flatten(run)
    if year == '2022pre':
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet' : "Summer22_22Sep2023_RunCD_V3_DATA_L1FastJet_AK4PFPuppi",
                'L2Relative' : "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi",
                'L3Absolute' : "Summer22_22Sep2023_RunCD_V3_DATA_L3Absolute_AK4PFPuppi",
                'L2L3Residual' : "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            # L2L3Residual Correction
            corr_L2L3 = evaluator[jec_names['L2L3Residual']].evaluate(run, eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3 * corr_L2L3
        ## MC Correction
        else:
            jec_names = {
                'L1FastJet' : "Summer22_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi",
                'L2Relative' : "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi",
                'L3Absolute' : "Summer22_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    elif year == '2022post':
        pass
    elif year == '2023pre':
        pass
    elif year == '2023post':
        pass
    if year == '2024':
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet' : "Summer24Prompt24_V2_DATA_L1FastJet_AK4PFPuppi",
                'L2Relative' : "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi",
                'L3Absolute' : "Summer24Prompt24_V2_DATA_L3Absolute_AK4PFPuppi",
                'L2L3Residual' : "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            # L2L3Residual Correction
            corr_L2L3 = evaluator[jec_names['L2L3Residual']].evaluate(run, eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3 * corr_L2L3
        ## MC Correction
        else:
            jec_names = {
                'L1FastJet' : "Summer24Prompt24_V2_MC_L1FastJet_AK4PFPuppi",
                'L2Relative' : "Summer24Prompt24_V2_MC_L2Relative_AK4PFPuppi",
                'L3Absolute' : "Summer24Prompt24_V2_MC_L3Absolute_AK4PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    return ak.unflatten(corr, counts)

def get_fjec_correction(year, pt, eta, phi, rho, area, run, isData):
    evaluator = correctionlib.CorrectionSet.from_file('data/JMESF/'+year+'/fatJet_jerc.json.gz')
    counts = ak.num(pt)
    run, _ = ak.broadcast_arrays(run, pt)
    rho, _ = ak.broadcast_arrays(rho, pt)
    pt, eta, phi, rho, area, run = ak.flatten(pt), ak.flatten(eta), ak.flatten(phi), ak.flatten(rho), ak.flatten(area), ak.flatten(run)
    if year == '2022pre':
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet' : "Summer22_22Sep2023_RunCD_V3_DATA_L1FastJet_AK8PFPuppi",
                'L2Relative' : "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK8PFPuppi",
                'L3Absolute' : "Summer22_22Sep2023_RunCD_V3_DATA_L3Absolute_AK8PFPuppi",
                'L2L3Residual' : "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK8PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            # L2L3Residual Correction
            corr_L2L3 = evaluator[jec_names['L2L3Residual']].evaluate(run, eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3 * corr_L2L3
        ## MC Correction
        else:
            jec_names = {
                'L1FastJet' : "Summer22_22Sep2023_V3_MC_L1FastJet_AK8PFPuppi",
                'L2Relative' : "Summer22_22Sep2023_V3_MC_L2Relative_AK8PFPuppi",
                'L3Absolute' : "Summer22_22Sep2023_V3_MC_L3Absolute_AK8PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    elif year == '2022post':
        pass
    elif year == '2023pre':
        pass
    elif year == '2023post':
        pass
    if year == '2024':
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet' : "Summer24Prompt24_V2_DATA_L1FastJet_AK8PFPuppi",
                'L2Relative' : "Summer24Prompt24_V2_DATA_L2Relative_AK8PFPuppi",
                'L3Absolute' : "Summer24Prompt24_V2_DATA_L3Absolute_AK8PFPuppi",
                'L2L3Residual' : "Summer24Prompt24_V2_DATA_L2L3Residual_AK8PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            # L2L3Residual Correction
            corr_L2L3 = evaluator[jec_names['L2L3Residual']].evaluate(run, eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3 * corr_L2L3
        ## MC Correction
        else:
            jec_names = {
                'L1FastJet' : "Summer24Prompt24_V2_MC_L1FastJet_AK8PFPuppi",
                'L2Relative' : "Summer24Prompt24_V2_MC_L2Relative_AK8PFPuppi",
                'L3Absolute' : "Summer24Prompt24_V2_MC_L3Absolute_AK8PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    return ak.unflatten(corr, counts)

####
# Muon ID scale factor
# https://twiki.cern.ch/twiki/bin/view/CMS/TWikiPAGsMUO
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonRun32022
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/MUO
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_mu_highpt_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_HighPtID_DEN_TrackerMuons"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

def get_mu_loose_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_LooseRelTkIso_DEN_HighPtID"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

def get_mu_hlt_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<52.),ak.full_like(pt,52.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

# All in one (Only three keys now)
def get_mu_sf (year, corr, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)

    if 'highpt' in corr:
        name = "NUM_HighPtID_DEN_TrackerMuons"
    elif 'iso' in corr: ## binning edge start 15 GeV
        name = "NUM_LooseRelTkIso_DEN_HighPtID"
    elif 'hlt' in corr: ## binning edge start 52 GeV
        name = "NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose"
        pt  = ak.where((pt<52.),ak.full_like(pt,52.),pt)
        flatpt = ak.flatten(pt)
    else:
        print("Wrong id")
    
    weight = evaluator[name].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

def get_mu_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

def get_mu_tight_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_TightID_DEN_TrackerMuons"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

def get_mu_tight_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_TightRelTkIso_DEN_HighPtID"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

###
# Muon scale and resolution (i.e. Rochester)
# https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
# RUN3 NOT UPDATED YET
###

####
# Photon ID scale factor
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_photon_id_sf(year, wp, eta, pt, phi):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/photon.json.gz')
    flateta, counts = ak.flatten(eta), ak.num(eta)
    pt  = ak.where((pt<20.),ak.full_like(pt,20.),pt)
    flatpt = ak.flatten(pt)
    flatphi = ak.flatten(phi)
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD',
        '2024': '2024'
    }
    if '2022' in year:
        sf_nominal = evaluator["Photon-ID-SF"].evaluate(yr[year], "sf", wp, flateta, flatpt)
        sf_up = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfup", wp, flateta, flatpt)
        sf_down = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfdown", wp, flateta, flatpt)
    elif '2023' in year:
        sf_nominal = evaluator["Photon-ID-SF"].evaluate(yr[year], "sf", wp, flateta, flatpt, flatphi)
        sf_up = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfup", wp, flateta, flatpt, flatphi)
        sf_down = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfdown", wp, flateta, flatpt, flatphi)
    elif '2024' in year:
        sf_nominal = evaluator["Photon-ID-SF"].evaluate(yr[year], "sf", wp, flateta, flatpt, flatphi)
        sf_up = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfup", wp, flateta, flatpt, flatphi)
        sf_down = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfdown", wp, flateta, flatpt, flatphi)
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)

####
# Electron ID scale factor
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_ele_veto_id_sf (year, eta, pt, phi):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD',
        '2024': '2024Prompt'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Veto", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Veto", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Veto", flateta, flatpt)
    elif '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Veto", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Veto", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Veto", flateta, flatpt, flatphi)
    elif '2024' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Veto", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Veto", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Veto", flateta, flatpt, flatphi)
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)

def get_ele_loose_id_sf (year, eta, pt, phi):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD',
        '2024': '2024Prompt'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Loose", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Loose", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Loose", flateta, flatpt)
    elif '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Loose", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Loose", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Loose", flateta, flatpt, flatphi)
    elif '2024' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Loose", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Loose", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Loose", flateta, flatpt, flatphi)
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)

def get_ele_medium_id_sf (year, eta, pt, phi):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Medium", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Medium", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Medium", flateta, flatpt)
    elif '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Medium", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Medium", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Medium", flateta, flatpt, flatphi)
    elif '2024' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Medium", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Medium", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Medium", flateta, flatpt, flatphi)
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)

def get_ele_tight_id_sf (year, eta, pt, phi):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD',
        '2024': '2024Prompt'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Tight", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Tight", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Tight", flateta, flatpt)
    elif '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Tight", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Tight", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Tight", flateta, flatpt, flatphi)
    elif '2024' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Tight", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Tight", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Tight", flateta, flatpt, flatphi)
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)

corrections = {}
corrections = {
    'get_pu_weight':            get_pu_weight,
    'get_met_xy_correction':    get_met_xy_correction,
    'get_jec_correction':       get_jec_correction,
    'get_fjec_correction':      get_fjec_correction,
    'get_mu_highpt_id_sf':      get_mu_highpt_id_sf,
    'get_mu_loose_iso_sf':      get_mu_loose_iso_sf,
    'get_mu_hlt_sf':            get_mu_hlt_sf,
    'get_mu_sf':                get_mu_sf,

    'get_mu_loose_id_sf':       get_mu_loose_id_sf,
    'get_mu_tight_id_sf':       get_mu_tight_id_sf,
    'get_mu_tight_iso_sf':      get_mu_tight_iso_sf,

    'get_photon_id_sf':         get_photon_id_sf,
    'get_ele_veto_id_sf':       get_ele_veto_id_sf,
    'get_ele_loose_id_sf':      get_ele_loose_id_sf,
    'get_ele_medium_id_sf':     get_ele_medium_id_sf,
    'get_ele_tight_id_sf':      get_ele_tight_id_sf,

#    'get_met_trig_weight':      get_met_trig_weight,
#    'get_ele_loose_id_sf':      get_ele_loose_id_sf,
#    'get_ele_tight_id_sf':      get_ele_tight_id_sf,
#    'get_ele_trig_weight':      get_ele_trig_weight,
#    'get_ele_reco_sf_below20':  get_ele_reco_sf_below20,
#    #'get_ele_reco_err_below20': get_ele_reco_err_below20,
#    'get_ele_reco_sf_above20':  get_ele_reco_sf_above20,
#    #'get_ele_reco_err_above20': get_ele_reco_err_above20,
#    'get_pho_tight_id_sf':      get_pho_tight_id_sf,
#    'get_pho_trig_weight':      get_pho_trig_weight,
#    'get_met_xy_correction':    XY_MET_Correction,
#    'get_nlo_ewk_weight':       get_nlo_ewk_weight,
#    'get_nnlo_nlo_weight':      get_nnlo_nlo_weight,
#    'get_ttbar_weight':         get_ttbar_weight,
#    'get_msd_corr':             get_msd_corr,
#    'get_btag_weight':          BTagCorrector,
#    'get_mu_rochester_sf':      get_mu_rochester_sf,
#    'jet_factory':              jet_factory,
#    'subjet_factory':           subjet_factory,
#    'fatjet_factory':           fatjet_factory,
#    'met_factory':              met_factory
}


save(corrections, 'data/corrections.coffea')
