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
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/LUM
####
#trueint = events.Pileup.nTrueInt
def get_pu_weight(year, trueint):
    correction = {
        '2022pre' : 'Collisions2022_355100_357900_eraBCD_GoldenJson',
        '2022post': 'Collisions2022_359022_362760_eraEFG_GoldenJson',
        '2023pre' : 'Collisions2023_366403_369802_eraBC_GoldenJson',
        '2023post': 'Collisions2023_369803_370790_eraD_GoldenJson',
    }
    evaluator = correctionlib.CorrectionSet.from_file('data/PUweight/'+year+'/puWeights.json.gz')
    weight = evaluator[correction[year]].evaluate(trueint, 'nominal')

    return weight

####
# MET XY Correction
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
####
def get_met_xy_correction(year, met_type, isData, met_pt, met_phi, npvGood):
    if year == '2022pre':
        epoch = '2022'
        yr = '2022'
    elif year == '2022post':
        epoch = '2022EE'
        yr = '2022'
    elif year == '2023pre':
        epoch = '2023'
        yr = '2023'
    elif year == '2023post':
        epoch = '2023BPix'
        yr = '2023'
    else:
        print("Now MET xy correction has only 2022, 2022EE, 2023, 2023BPix")
    evaluator = correctionlib.CorrectionSet.from_file('data/JetMETCorr/'+year+'/met_xyCorrections_'+yr+'_'+epoch+'.json.gz')

    if isData:
        dtmc = 'DATA'
    else:
        dtmc = 'MC'

    corr_pt  = evaluator['met_xy_corrections'].evaluate('pt',  met_type, epoch, dtmc, 'nom', met_pt, met_phi, npvGood)
    corr_phi = evaluator['met_xy_corrections'].evaluate('phi', met_type, epoch, dtmc, 'nom', met_pt, met_phi, npvGood)

    return corr_pt, corr_phi

####
# JEC
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
####
def get_jec_correction(year, pt, eta, phi, rho, area, run, isData):
    jecTag = {
        "2022pre":  "Summer22_22Sep2023",
        "2022post": "Summer22EE_22Sep2023",
        "2023pre":  "Summer23Prompt23",
        "2023post": "Summer23BPixPrompt23",
        }
    evaluator = correctionlib.CorrectionSet.from_file('data/JetMETCorr/'+year+'/jet_jerc.json.gz')
    counts = ak.num(pt)
    run, _ = ak.broadcast_arrays(run, pt)
    rho, _ = ak.broadcast_arrays(rho, pt)
    pt, eta, phi, rho, area, run = ak.flatten(pt), ak.flatten(eta), ak.flatten(phi), ak.flatten(rho), ak.flatten(area), ak.flatten(run)
    if '2022' in year:
        ## DATA Correction
        if isData:
            if run[0] >= 355794 and run[0] <= 359021:
                dataArea = "RunCD"
            if run[0] >= 359022 and run[0] <= 360331:
                dataArea = "RunE"
            if run[0] >= 360332 and run[0] <= 362180:
                dataArea = "RunF"
            if run[0] >= 362350 and run[0] <= 362760:
                dataArea = "RunG"

            jec_names = {
                'L1FastJet'   : f"{jecTag[year]}_V3_DATA_L1FastJet_AK4PFPuppi",
                'L2Relative'  : f"{jecTag[year]}_V3_DATA_L2Relative_AK4PFPuppi",
                'L3Absolute'  : f"{jecTag[year]}_V3_DATA_L3Absolute_AK4PFPuppi",
                'L2L3Residual': f"{jecTag[year]}_V3_DATA_L2L3Residual_AK4PFPuppi"
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
                'L1FastJet'  : f"{jecTag[year]}_V3_MC_L1FastJet_AK4PFPuppi",
                'L2Relative' : f"{jecTag[year]}_V3_MC_L2Relative_AK4PFPuppi",
                'L3Absolute' : f"{jecTag[year]}_V3_MC_L3Absolute_AK4PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    if '2023' in year:
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet'   : f"{jecTag[year]}_V3_DATA_L1FastJet_AK4PFPuppi",
                'L2Relative'  : f"{jecTag[year]}_V3_DATA_L2Relative_AK4PFPuppi",
                'L3Absolute'  : f"{jecTag[year]}_V3_DATA_L3Absolute_AK4PFPuppi",
                'L2L3Residual': f"{jecTag[year]}_V3_DATA_L2L3Residual_AK4PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            if 'post' in year:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            else:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            # L2L3Residual Correction
            corr_L2L3 = evaluator[jec_names['L2L3Residual']].evaluate(run, eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3 * corr_L2L3
        ## MC Correction
        else:
            jec_names = {
                'L1FastJet'  : f"{jecTag[year]}_V3_MC_L1FastJet_AK4PFPuppi",
                'L2Relative' : f"{jecTag[year]}_V3_MC_L2Relative_AK4PFPuppi",
                'L3Absolute' : f"{jecTag[year]}_V3_MC_L3Absolute_AK4PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            if 'post' in year:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            else:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    if year == '2024':
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet'   : "Summer24Prompt24_V2_DATA_L1FastJet_AK4PFPuppi",
                'L2Relative'  : "Summer24Prompt24_V2_DATA_L2Relative_AK4PFPuppi",
                'L3Absolute'  : "Summer24Prompt24_V2_DATA_L3Absolute_AK4PFPuppi",
                'L2L3Residual': "Summer24Prompt24_V2_DATA_L2L3Residual_AK4PFPuppi"
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
    jecTag = {
        "2022pre":  "Summer22_22Sep2023",
        "2022post": "Summer22EE_22Sep2023",
        "2023pre":  "Summer23Prompt23",
        "2023post": "Summer23BPixPrompt23",
        }
    evaluator = correctionlib.CorrectionSet.from_file('data/JetMETCorr/'+year+'/fatJet_jerc.json.gz')
    counts = ak.num(pt)
    run, _ = ak.broadcast_arrays(run, pt)
    rho, _ = ak.broadcast_arrays(rho, pt)
    pt, eta, phi, rho, area, run = ak.flatten(pt), ak.flatten(eta), ak.flatten(phi), ak.flatten(rho), ak.flatten(area), ak.flatten(run)
    
    if '2022' in year:
        ## DATA Correction
        if isData:
            if run[0] >= 355794 and run[0] <= 359021:
                dataArea = "RunCD"
            if run[0] >= 359022 and run[0] <= 360331:
                dataArea = "RunE"
            if run[0] >= 360332 and run[0] <= 362180:
                dataArea = "RunF"
            if run[0] >= 362350 and run[0] <= 362760:
                dataArea = "RunG"

            jec_names = {
                'L1FastJet'   : f"{jecTag[year]}_V3_DATA_L1FastJet_AK8PFPuppi",
                'L2Relative'  : f"{jecTag[year]}_V3_DATA_L2Relative_AK8PFPuppi",
                'L3Absolute'  : f"{jecTag[year]}_V3_DATA_L3Absolute_AK8PFPuppi",
                'L2L3Residual': f"{jecTag[year]}_V3_DATA_L2L3Residual_AK8PFPuppi"
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
                'L1FastJet'  : f"{jecTag[year]}_V3_MC_L1FastJet_AK8PFPuppi",
                'L2Relative' : f"{jecTag[year]}_V3_MC_L2Relative_AK8PFPuppi",
                'L3Absolute' : f"{jecTag[year]}_V3_MC_L3Absolute_AK8PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta,  pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    if '2023' in year:
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet'   : f"{jecTag[year]}_V3_DATA_L1FastJet_AK8PFPuppi",
                'L2Relative'  : f"{jecTag[year]}_V3_DATA_L2Relative_AK8PFPuppi",
                'L3Absolute'  : f"{jecTag[year]}_V3_DATA_L3Absolute_AK8PFPuppi",
                'L2L3Residual': f"{jecTag[year]}_V3_DATA_L2L3Residual_AK8PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            if 'post' in year:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            else:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            # L2L3Residual Correction
            corr_L2L3 = evaluator[jec_names['L2L3Residual']].evaluate(run, eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3 * corr_L2L3
        ## MC Correction
        else:
            jec_names = {
                'L1FastJet'  : f"{jecTag[year]}_V3_MC_L1FastJet_AK8PFPuppi",
                'L2Relative' : f"{jecTag[year]}_V3_MC_L2Relative_AK8PFPuppi",
                'L3Absolute' : f"{jecTag[year]}_V3_MC_L3Absolute_AK8PFPuppi"
            }
            # L1FastJet Correction
            corr_L1 = evaluator[jec_names['L1FastJet']].evaluate(area, eta, pt, rho)
            # L2Relative Correction
            if 'post' in year:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, phi, pt)
            else:
                corr_L2 = evaluator[jec_names['L2Relative']].evaluate(eta, pt)
            # L3Absolute Correction
            corr_L3 = evaluator[jec_names['L3Absolute']].evaluate(eta, pt)
            corr = corr_L1 * corr_L2 * corr_L3

    if year == '2024':
        ## DATA Correction
        if isData:
            jec_names = {
                'L1FastJet'   : "Summer24Prompt24_V2_DATA_L1FastJet_AK8PFPuppi",
                'L2Relative'  : "Summer24Prompt24_V2_DATA_L2Relative_AK8PFPuppi",
                'L3Absolute'  : "Summer24Prompt24_V2_DATA_L3Absolute_AK8PFPuppi",
                'L2L3Residual': "Summer24Prompt24_V2_DATA_L2L3Residual_AK8PFPuppi"
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
                'L1FastJet'  : "Summer24Prompt24_V2_MC_L1FastJet_AK8PFPuppi",
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

def get_mu_loose_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_LoosePFIso_DEN_LooseID"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

def get_mu_tight_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_TightPFIso_DEN_TightID"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

# High pT muon (For Nanoaod Study)
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
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Photon-ID-SF"].evaluate(yr[year], "sf",     wp, flateta, flatpt)
        sf_up      = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfup",   wp, flateta, flatpt)
        sf_down    = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfdown", wp, flateta, flatpt)
    if '2023' in year:
        sf_nominal = evaluator["Photon-ID-SF"].evaluate(yr[year], "sf",     wp, flateta, flatpt, flatphi)
        sf_up      = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfup",   wp, flateta, flatpt, flatphi)
        sf_down    = evaluator["Photon-ID-SF"].evaluate(yr[year], "sfdown", wp, flateta, flatpt, flatphi)
    
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
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Veto", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Veto", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Veto", flateta, flatpt)
    if '2023' in year:
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
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Loose", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Loose", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Loose", flateta, flatpt)
    if '2023' in year:
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
    if '2023' in year:
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
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Tight", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Tight", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Tight", flateta, flatpt)
    if '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Tight", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Tight", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Tight", flateta, flatpt, flatphi)
    
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)


####
# Electron Reco scale factor
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammSFandSSRun3
####

def get_ele_reco_sf_below20(year, eta, pt, phi):

    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    pt = ak.where((pt>19.99), ak.full_like(pt,19.99), pt)

    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "RecoBelow20", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "RecoBelow20", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "RecoBelow20", flateta, flatpt)
    if '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "RecoBelow20", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "RecoBelow20", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "RecoBelow20", flateta, flatpt, flatphi)
    
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)


def get_ele_reco_sf_20to75(year, eta, pt, phi):

    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')
    
    pt = ak.where((pt<20.), ak.full_like(pt,20.), pt)
    pt = ak.where((pt>74.99), ak.full_like(pt,74.99), pt)

    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Reco20to75", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Reco20to75", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Reco20to75", flateta, flatpt)
    if '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Reco20to75", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "Reco20to75", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "Reco20to75", flateta, flatpt, flatphi)
    
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)


def get_ele_reco_sf_Above75(year, eta, pt, phi):

    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')
    
    pt = ak.where((pt<75.), ak.full_like(pt,75.), pt)

    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "RecoAbove75", flateta, flatpt)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "RecoAbove75", flateta, flatpt)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "RecoAbove75", flateta, flatpt)
    if '2023' in year:
        sf_nominal = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "RecoAbove75", flateta, flatpt, flatphi)
        sf_up = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfup", "RecoAbove75", flateta, flatpt, flatphi)
        sf_down = evaluator["Electron-ID-SF"].evaluate(yr[year], "sfdown", "RecoAbove75", flateta, flatpt, flatphi)
    
    return ak.unflatten(sf_nominal, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)


###
# MET trigger efficiency SFs. Depends on recoil.
###
#
def get_met_trig_weight(year, met):
    corrname = {
        '2022pre' : "recoil_trigger_sf_22pre",
        '2022post': "recoil_trigger_sf_22post",
        '2023pre' : "recoil_trigger_sf_23pre",
        '2023post': "recoil_trigger_sf_23post"
    }
    cset = correctionlib.CorrectionSet.from_file(f"data/METTrigEff/recoil_trigger_sf_{year}.json.gz")
    corr = cset[corrname[year]]

    met  = ak.fill_none(met, 0.)

    sf = corr.evaluate(met, "nominal")
    sf_up   = corr.evaluate(met, "up")
    sf_down = corr.evaluate(met, "down")

    return sf, sf_up, sf_down



####
# Electron Trigger weight
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3#Scale_factors_and_correction_AN2
# edges eta [-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5]
#        pt [25.0, 30.0, 32.0, 35.0, 40.0, 45.0, 50.0, 75.0, 100.0, 200.0, 500.0]
####

def get_ele_trig_weight(year, eta, pt, cutbased):
    evaluator = correctionlib.CorrectionSet.from_file("data/EGammaSF/"+year+"/electronHlt.json.gz")
    yr = {
    '2022pre' : '2022Re-recoBCD',
    '2022post': '2022Re-recoE+PromptFG',
    '2023pre' : '2023PromptC',
    '2023post': '2023PromptD'
    }

    wp = "HLT_SF_Ele30_"+str(cutbased)+"ID"

    pt = ak.where((pt<25.0), ak.full_like(pt,25.), pt)

    flatpt = ak.flatten(pt)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    sf = evaluator["Electron-HLT-SF"].evaluate(yr[year], "sf", wp, flateta, flatpt)
    sf_up   = evaluator["Electron-HLT-SF"].evaluate(yr[year], "sfup",   wp, flateta, flatpt)
    sf_down = evaluator["Electron-HLT-SF"].evaluate(yr[year], "sfdown", wp, flateta, flatpt)

    return ak.unflatten(sf, counts=counts), ak.unflatten(sf_up, counts=counts), ak.unflatten(sf_down, counts=counts)

####
# Photon Trigger weight
# Measured using MET, GJ (trg: Photon200, ref: PFHT1050)
####

def get_pho_trig_weight(year, pt):
    corrname = {
        "2022pre" : "photon_trigger_sf_2022pre",
        "2022post": "photon_trigger_sf_2022post",
        "2023pre" : "photon_trigger_sf_2023pre",
        "2023post": "photon_trigger_sf_2023post",
    }

    cset = correctionlib.CorrectionSet.from_file(f"data/PhotonTrigEff/photon_trigger_sf_{year}.json.gz")
    corr = corr = cset[corrname[year]]

    pt = ak.fill_none(pt, 0.)

    sf      = corr.evaluate(pt, "nominal")
    sf_up   = corr.evaluate(pt, "up")
    sf_down = corr.evaluate(pt, "down")


    return sf, sf_up, sf_down



def get_ttbar_weight(pt):
   return (0.103*np.exp(-0.0118 * np.clip(pt, 0, 800))-0.000134*np.clip(pt, 0, 800)+0.973)*(0.991+0.000075*np.clip(pt, 0, 800))

###
# V+jets NLO k-factors
# Combination of NNLO QCD and NLO EWK
###
def get_nnlo_nlo_wjet_AN24_075(channel, mass):
    # The W background K factors and uncertainties for electron and muon channel.
    # Combination of NNLO QCD and NLO EWK is done with additive + mixed term approach.
    # Provided by AN2024_075_v11.
    kfactors = {
        "electron": [
            (120, 200, 1.143),
            (200, 400, 1.217),
            (400, 800, 1.215),
            (800, 1500, 1.214),
            (1500, 2500, 1.168),
            (2500, 4000, 1.148),
            (4000, 6000, 1.102),
            (6000, 8000, 1.084),
        ],
        "muon": [
            (120, 200, 1.112),
            (200, 400, 1.165),
            (400, 800, 1.161),
            (800, 1500, 1.152),
            (1500, 2500, 1.100),
            (2500, 4000, 1.084),
            (4000, 6000, 1.050),
            (6000, 8000, 1.040),
        ],
    }

    edges = np.array([low for low, _, _ in kfactors[channel]] + [8000])
    values = np.array([v for _, _, v in kfactors[channel]])

    m = ak.to_numpy(ak.flatten(mass))
    m = np.clip(m, edges[0], edges[-1] - 1e-6)

    k = np.interp(m, edges[:-1], values)
    return ak.unflatten(k, ak.num(mass, axis=-1))


###
# V+jets NLO k-factors
# Only use nlo ewk sf
###

nlo_ewk_hists = {
        'dy': ["* * data/Run2Files/vjets_SFs/merged_kfactors_zjets.root"],
        'w':  ["* * data/Run2Files/vjets_SFs/merged_kfactors_wjets.root"],
        'z':  ["* * data/Run2Files/vjets_SFs/merged_kfactors_zjets.root"],
        'a':  ["* * data/Run2Files/vjets_SFs/merged_kfactors_gjets.root"],
}    
get_nlo_ewk_weight = {}
for p in ['dy','w','z','a']:
    print(nlo_ewk_hists[p])
    ext = extractor()
    ext.add_weight_sets(nlo_ewk_hists[p])
    ext.finalize()
    get_nlo_ewk_weight[p] = ext.make_evaluator()["kfactor_monojet_ewk"]

####
## V+jets NNLO weights
## The schema is process_NNLO_NLO_QCD1QCD2QCD3_EW1EW2EW3_MIX, where 'n' stands for 'nominal', 'u' for 'up', and 'd' for 'down'
####
#
#histname={
#    'dy': 'eej_NNLO_NLO_',
#    'w':  'evj_NNLO_NLO_',
#    'z': 'vvj_NNLO_NLO_',
#    'a': 'aj_NNLO_NLO_'
#}
#correlated_variations = {
#    'cen':    'nnn_nnn_n',
#    'qcd1up': 'unn_nnn_n',
#    'qcd1do': 'dnn_nnn_n',
#    'qcd2up': 'nun_nnn_n',
#    'qcd2do': 'ndn_nnn_n',
#    'qcd3up': 'nnu_nnn_n',
#    'qcd3do': 'nnd_nnn_n',
#    'ew1up' : 'nnn_unn_n',
#    'ew1do' : 'nnn_dnn_n',
#    'mixup' : 'nnn_nnn_u',
#    'mixdo' : 'nnn_nnn_d',
#    'muFup' : 'nnn_nnn_n_Weight_scale_variation_muR_1p0_muF_2p0',
#    'muFdo' : 'nnn_nnn_n_Weight_scale_variation_muR_1p0_muF_0p5',
#    'muRup' : 'nnn_nnn_n_Weight_scale_variation_muR_2p0_muF_1p0',
#    'muRdo' : 'nnn_nnn_n_Weight_scale_variation_muR_0p5_muF_1p0'
#}
#uncorrelated_variations = {
#    'dy': {
#        'ew2Gup': 'nnn_nnn_n',
#        'ew2Gdo': 'nnn_nnn_n',
#        'ew2Wup': 'nnn_nnn_n',
#        'ew2Wdo': 'nnn_nnn_n',
#        'ew2Zup': 'nnn_nun_n',
#        'ew2Zdo': 'nnn_ndn_n',
#        'ew3Gup': 'nnn_nnn_n',
#        'ew3Gdo': 'nnn_nnn_n',
#        'ew3Wup': 'nnn_nnn_n',
#        'ew3Wdo': 'nnn_nnn_n',
#        'ew3Zup': 'nnn_nnu_n',
#        'ew3Zdo': 'nnn_nnd_n'
#    },
#    'w': {
#        'ew2Gup': 'nnn_nnn_n',
#        'ew2Gdo': 'nnn_nnn_n',
#        'ew2Wup': 'nnn_nun_n',
#        'ew2Wdo': 'nnn_ndn_n',
#        'ew2Zup': 'nnn_nnn_n',
#        'ew2Zdo': 'nnn_nnn_n',
#        'ew3Gup': 'nnn_nnn_n',
#        'ew3Gdo': 'nnn_nnn_n',
#        'ew3Wup': 'nnn_nnu_n',
#        'ew3Wdo': 'nnn_nnd_n',
#        'ew3Zup': 'nnn_nnn_n',
#        'ew3Zdo': 'nnn_nnn_n'
#    },
#    'z': {
#        'ew2Gup': 'nnn_nnn_n',
#        'ew2Gdo': 'nnn_nnn_n',
#        'ew2Wup': 'nnn_nnn_n',
#        'ew2Wdo': 'nnn_nnn_n',
#        'ew2Zup': 'nnn_nun_n',
#        'ew2Zdo': 'nnn_ndn_n',
#        'ew3Gup': 'nnn_nnn_n',
#        'ew3Gdo': 'nnn_nnn_n',
#        'ew3Wup': 'nnn_nnn_n',
#        'ew3Wdo': 'nnn_nnn_n',
#        'ew3Zup': 'nnn_nnu_n',
#        'ew3Zdo': 'nnn_nnd_n'
#    },
#    'a': {
#        'ew2Gup': 'nnn_nun_n',
#        'ew2Gdo': 'nnn_ndn_n',
#        'ew2Wup': 'nnn_nnn_n',
#        'ew2Wdo': 'nnn_nnn_n',
#        'ew2Zup': 'nnn_nnn_n',
#        'ew2Zdo': 'nnn_nnn_n',
#        'ew3Gup': 'nnn_nnu_n',
#        'ew3Gdo': 'nnn_nnd_n',
#        'ew3Wup': 'nnn_nnn_n',
#        'ew3Wdo': 'nnn_nnn_n',
#        'ew3Zup': 'nnn_nnn_n',
#        'ew3Zdo': 'nnn_nnn_n'
#    }
#}
#get_nnlo_nlo_weight = {}
#for year in ['2016postVFP', '2016preVFP', '2017','2018']:
#    get_nnlo_nlo_weight[year] = {}
#    if '2016' in year:
#        nnlo_file = {
#            'dy': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_eej_madgraph_2016.root"),
#            'w': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_evj_madgraph_2016.root"),
#            'z': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_vvj_madgraph_2016.root"),
#            'a': uproot.open("data/Vboson_Pt_Reweighting/2016/TheoryXS_aj_madgraph_2016.root")
#        }
#    else:
#        nlo_file = {
#            'dy': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_eej_madgraph_"+year+".root"),
#            'w': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_evj_madgraph_"+year+".root"),
#            'z': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_vvj_madgraph_"+year+".root"),
#            'a': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_aj_madgraph_"+year+".root")
#        }
#    for p in ['dy','w','z','a']:
#        get_nnlo_nlo_weight[year][p] = {}
#        for cv in correlated_variations:
#            histo=nnlo_file[p][histname[p]+correlated_variations[cv]]
#            get_nnlo_nlo_weight[year][p][cv]=lookup_tools.dense_lookup.dense_lookup(histo.values(), histo.axes)
#        for uv in uncorrelated_variations[p]:
#            histo=nnlo_file[p][histname[p]+uncorrelated_variations[p][uv]]
#            get_nnlo_nlo_weight[year][p][uv]=lookup_tools.dense_lookup.dense_lookup(histo.values(), histo.axes)
#
#
#
#
## Soft drop mass correction updated for UL. Copied from:
## https://github.com/jennetd/hbb-coffea/blob/master/boostedhiggs/corrections.py
## Renamed copy of corrected_msoftdrop() function
#
#msdcorr = correctionlib.CorrectionSet.from_file('data/msdcorr.json')
#
#def get_msd_corr(fatjets):
#    msdraw = np.sqrt(
#        np.maximum(
#            0.0,
#            (fatjets.subjets * (1 - fatjets.subjets.rawFactor)).sum().mass2,
#        )
#    )
#    msdfjcorr = msdraw / (1 - fatjets.rawFactor)
#
#    corr = msdcorr["msdfjcorr"].evaluate(
#        np.array(ak.flatten(msdfjcorr / fatjets.pt)),
#        np.array(ak.flatten(np.log(fatjets.pt))),
#        np.array(ak.flatten(fatjets.eta)),
#    )
#    corr = ak.unflatten(corr, ak.num(fatjets))
#    corrected_mass = msdfjcorr * corr
#
#    return corrected_mass
#

from coffea.lookup_tools.correctionlib_wrapper import correctionlib_wrapper
from coffea.lookup_tools.dense_lookup import dense_lookup

class BTagCorrector:

    def __init__(self, tagger, year, workingpoint, caller):
        self._year = year

        wp = {}
        wp['loose'] = 'L'
        wp['medium'] = 'M'
        wp['tight'] = 'T'
        wp['verytight'] = 'XT'
        wp['veryverytight'] = 'XXT'
        self._wp = wp[workingpoint]
        self._mc = caller 

        btvjson = {}
        if year == 2024:
            btvjson['PNetUParT'] = {
                'comb': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["UParTAK4_comb"],
                'ligh': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["UParTAK4_light"],
            }
        else:
            btvjson['PNetUParT'] = {
                'comb': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["particleNet_comb"],
                'ligh': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["particleNet_light"],
            }
        self.sf = btvjson[tagger]

        files = {
            '2022pre' : 'btageff2022pre.merged',
            '2022post': 'btageff2022post.merged',
            '2023pre' : 'btageff2023pre.merged',
            '2023post': 'btageff2023post.merged',
        }
        filename = 'hists/'+files[year]
        btag_file = load(filename)
        #for k in btag_file[tagger]:
        #    try:
        #        btag += btag_file[tagger][k]
        #    except:
        #        btag = btag_file[tagger][k]
        btag = btag_file[tagger][caller]
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
        pt = ak.where((pt<20.0), ak.full_like(pt,20.0), pt)
        flatpt =  ak.fill_none(ak.flatten(pt), 20.)

        flatflavor = ak.fill_none(ak.flatten(flavor), 0)
        
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
        def P(eff):
            weight = ak.where(istag, eff, 1-eff)
            return ak.prod(weight, axis=1)

        '''
        particleNet_comb
         =>  systematic ,
         =>  working_point ,  L/M/T/XT/XXT
         =>  flavor ,  hadron flavor definition: 5=b, 4=c, 0=udsg
         =>  abseta ,
         =>  pt ,
        '''

        eff = ak.where(
            ~np.isnan(ak.fill_none(pt, np.nan)),
            ak.unflatten(self.eff.evaluate(flatflavor, flatpt, flateta), counts=counts),
            ak.zeros_like(pt)
        )
        sf_nom = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_up_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_down_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_up_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_down_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('down_uncorrelated',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('down_uncorrelated',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)        
            )
        )
        sf_light_up_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_down_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_up_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_down_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('down_uncorrelated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
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

###
# Muon scale and resolution (i.e. Rochester)
# https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
# RUN3 NOT UPDATED YET
###
#
#tag = 'roccor.Run2.v5'
#get_mu_rochester_sf = {}
#for year in ['2016postVFP', '2016preVFP', '2017','2018']:
#    if '2016postVFP' in year: 
#        fname = f'data/{tag}/RoccoR2016bUL.txt'
#    elif '2016preVFP' in year:  
#        fname = f'data/{tag}/RoccoR2016aUL.txt'
#    else:
#        fname = f'data/{tag}/RoccoR{year}UL.txt'
#    sfs = lookup_tools.txt_converters.convert_rochester_file(fname,loaduncs=True)
#    get_mu_rochester_sf[year] = lookup_tools.rochester_lookup.rochester_lookup(sfs)
#



corrections = {}
corrections = {
    'get_pu_weight':            get_pu_weight,
    'get_met_xy_correction':    get_met_xy_correction,
    'get_jec_correction':       get_jec_correction,
    'get_fjec_correction':      get_fjec_correction,

    'get_mu_sf':                get_mu_sf,

    'get_mu_loose_id_sf':       get_mu_loose_id_sf,
    'get_mu_tight_id_sf':       get_mu_tight_id_sf,
    'get_mu_loose_iso_sf':      get_mu_loose_iso_sf,
    'get_mu_tight_iso_sf':      get_mu_tight_iso_sf,

    'get_photon_id_sf':         get_photon_id_sf,
    'get_ele_loose_id_sf':      get_ele_loose_id_sf,
    'get_ele_tight_id_sf':      get_ele_tight_id_sf,

    'get_ele_reco_sf_below20':  get_ele_reco_sf_below20,
    'get_ele_reco_sf_20to75':  get_ele_reco_sf_20to75,
    'get_ele_reco_sf_Above75':  get_ele_reco_sf_Above75,

    'get_met_trig_weight':      get_met_trig_weight,
    'get_ele_trig_weight':      get_ele_trig_weight,
    'get_pho_trig_weight':      get_pho_trig_weight,
    
    'get_ttbar_weight':         get_ttbar_weight,

    'get_btag_weight':          BTagCorrector,

    'get_nlo_ewk_weight':       get_nlo_ewk_weight,
#    'get_nnlo_nlo_wjet':         get_nnlo_nlo_wjet,
#    'get_ele_medium_id_sf':     get_ele_medium_id_sf,
#    'get_ele_veto_id_sf':       get_ele_veto_id_sf,
#    'get_nnlo_nlo_weight':      get_nnlo_nlo_weight,
#    'get_msd_corr':             get_msd_corr,
#    'get_mu_rochester_sf':      get_mu_rochester_sf,
}


save(corrections, 'data/corrections.coffea')
