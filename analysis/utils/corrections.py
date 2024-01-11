#! /usr/bin/env python
import correctionlib
import os
#import uproot
#import uproot_methods
import awkward as ak

import numpy as np
from coffea import hist, lookup_tools
from coffea.lookup_tools import extractor, dense_lookup

import uproot3
from coffea import util
from coffea.util import save, load
import json
###
# MET trigger efficiency SFs, 2017/18 from monojet. Depends on recoil.
###

met_trig_hists = {
    '2016': uproot3.open("data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")['hden_monojet_recoil_clone_passed'],
    '2017': uproot3.open("data/trigger_eff/met_trigger_sf.root")['120pfht_hltmu_1m_2017'],
    '2018': uproot3.open("data/trigger_eff/met_trigger_sf.root")['120pfht_hltmu_1m_2018']
}
get_met_trig_weight = {}
for year in ['2016','2017','2018']:
    met_trig_hist=met_trig_hists[year]
    get_met_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(met_trig_hist.values, met_trig_hist.edges)


####
# Electron ID scale factor
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_ele_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    flatpt = ak.flatten(pt)
    weight = evaluator["UL-Electron-ID-SF"].evaluate(year, "sf", "Loose", flateta, flatpt)

    return ak.unflatten(weight, counts=counts)

def get_ele_tight_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    flatpt = ak.flatten(pt)
    weight = evaluator["UL-Electron-ID-SF"].evaluate(year, "sf", "Tight", flateta, flatpt)

    return ak.unflatten(weight, counts=counts)


####
# Electron Trigger weight
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
# Copy from previous correctionsUL.py file
####

ele_trig_hists = {
    '2016postVFP': uproot3.open("data/ElectronTrigEff/egammaEffi.txt_EGM2D-2016postVFP.root")['EGamma_SF2D'],
    '2016preVFP' : uproot3.open("data/ElectronTrigEff/egammaEffi.txt_EGM2D-2016preVFP.root")['EGamma_SF2D'],
    '2017': uproot3.open("data/ElectronTrigEff/egammaEffi.txt_EGM2D-2017.root")['EGamma_SF2D'],#monojet measurement for the combined trigger path
    '2018': uproot3.open("data/ElectronTrigEff/egammaEffi.txt_EGM2D-2018.root")['EGamma_SF2D'] #approved by egamma group: https://indico.cern.ch/event/924522/
}
get_ele_trig_weight = {}
for year in ['2016postVFP', '2016preVFP', '2017','2018']:
    ele_trig_hist = ele_trig_hists[year]
    get_ele_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(ele_trig_hist.values, ele_trig_hist.edges)



####
# Electron Reco scale factor
# root files: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
# Code Copy from previous correctionsUL.py file
####

ele_reco_files_below20 = {
    '2016postVFP': uproot3.open("data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root"),
    '2016preVFP': uproot3.open("data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root"),
    '2017': uproot3.open("data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2017.root"),
    '2018': uproot3.open("data/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root")
}
get_ele_reco_sf_below20 = {}
get_ele_reco_err_below20 = {}
for year in ['2016postVFP','2016preVFP','2017','2018']:
    ele_reco_hist = ele_reco_files_below20[year]["EGamma_SF2D"]
    get_ele_reco_sf_below20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.values, ele_reco_hist.edges)
    get_ele_reco_err_below20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.variances ** 0.5, ele_reco_hist.edges)


ele_reco_files_above20 = {
    '2016postVFP': uproot3.open("data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root"),
    '2016preVFP': uproot3.open("data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root"),
    '2017': uproot3.open("data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root"),
    '2018': uproot3.open("data/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root")
}
get_ele_reco_sf_above20 = {}
get_ele_reco_err_above20 = {}
for year in ['2016postVFP','2016preVFP','2017','2018']:
    ele_reco_hist = ele_reco_files_above20[year]["EGamma_SF2D"]
    get_ele_reco_sf_above20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.values, ele_reco_hist.edges)
    get_ele_reco_err_above20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.variances ** 0.05, ele_reco_hist.edges)
    


####
# Photon ID scale factor
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_pho_tight_id_sf(year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/photon.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    flatpt = ak.flatten(pt)
    weight = evaluator["UL-Photon-ID-SF"].evaluate(year, "sf", "Tight", flateta, flatpt)

    return ak.unflatten(weight, counts=counts)


def get_pho_loose_id_sf(year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/photon.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
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

pho_trig_files = {
    '2016postVFP': uproot3.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"),
    '2016preVFP': uproot3.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"),
    "2017": uproot3.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"),
    "2018": uproot3.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
}
get_pho_trig_weight = {}
for year in ['2016postVFP','2016preVFP','2017','2018']:
    pho_trig_hist = pho_trig_files[year]["hden_photonpt_clone_passed"]
    get_pho_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(pho_trig_hist.values, pho_trig_hist.edges)



# https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL

####
# Muon ID scale factor
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018n?topic=MuonUL2018
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/MUO
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_muon_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/muon_Z.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    flatpt = ak.flatten(pt)
    if year == '2018':
        weight = evaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(year, flateta, flatpt, "sf")
    else:
        weight = evaluator["NUM_LooseID_DEN_genTracks"].evaluate(year, flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)

def get_muon_tight_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/muon_Z.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    flatpt = ak.flatten(pt)
    if year == '2018':
        weight = evaluator["NUM_TightID_DEN_TrackerMuons"].evaluate(year, flateta, flatpt, "sf")
    else:
        weight = evaluator["NUM_TightID_DEN_genTracks"].evaluate(year, flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)




####
# Muon Iso scale factor
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018n?topic=MuonUL2018
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/MUO
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_muon_loose_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/muon_Z.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    flatpt = ak.flatten(pt)
    weight = evaluator["NUM_LooseRelIso_DEN_LooseID"].evaluate(year, flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)

def get_muon_tight_iso_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'_UL/muon_Z.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    flatpt = ak.flatten(pt)
    weight = evaluator["NUM_TightRelIso_DEN_TightIDandIPCut"].evaluate(year, flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)


####
# PU weight
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/LUM
####
#trueint = events.Pileup.nTrueInt
def pu_weight(year, trueint):
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

def XY_MET_Correction(year, events, pt, phi):
    if 'genWeight' in events.fields:
        isMC = True
    else:
        isData = True

    npv = events.PV.npvsGood
    mask = np.asarray(npv>100)
    npv = np.asarray(npv)
    npv[mask==True] = 100
    run = events.run

    evaluator = correctionlib.CorrectionSet.from_file('data/JetMETCorr/'+year+'_UL/met.json.gz')

    if isData:
        corrected_pt = evaluator['pt_metphicorr_pfmet_data'].evaluate(pt,phi,npv,run)
        corrected_phi = evaluator['phi_metphicorr_pfmet_data'].evaluate(pt,phi,npv,run)

    if isMC:
        corrected_pt = evaluator['pt_metphicorr_pfmet_mc'].evaluate(pt,phi,npv,run)
        corrected_phi = evaluator['phi_metphicorr_pfmet_mc'].evaluate(pt,phi,npv,run)

    return corrected_pt, corrected_phi



###
# V+jets NLO k-factors
# Only use nlo ewk sf
###

nlo_ewk_hists = {
    'dy': uproot3.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_ewk"],
    'w': uproot3.open("data/vjets_SFs/merged_kfactors_wjets.root")["kfactor_monojet_ewk"],
    'z': uproot3.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_ewk"],
    'a': uproot3.open("data/vjets_SFs/merged_kfactors_gjets.root")["kfactor_monojet_ewk"]
}    
get_nlo_ewk_weight = {}
for year in ['2016','2017','2018']:
    get_nlo_ewk_weight[year] = {}
    for p in ['dy','w','z','a']:
        get_nlo_ewk_weight[year][p] = lookup_tools.dense_lookup.dense_lookup(nlo_ewk_hists[p].values, nlo_ewk_hists[p].edges)



def get_ttbar_weight(pt):
    return np.exp(0.0615 - 0.0005 * np.clip(pt, 0, 800))

def get_msd_weight(pt, eta):
    gpar = np.array([1.00626, -1.06161, 0.0799900, 1.20454])
    cpar = np.array([1.09302, -0.000150068, 3.44866e-07, -2.68100e-10, 8.67440e-14, -1.00114e-17])
    fpar = np.array([1.27212, -0.000571640, 8.37289e-07, -5.20433e-10, 1.45375e-13, -1.50389e-17])
    genw = gpar[0] + gpar[1]*np.power(pt*gpar[2], -gpar[3])
    ptpow = np.power.outer(pt, np.arange(cpar.size))
    cenweight = np.dot(ptpow, cpar)
    forweight = np.dot(ptpow, fpar)
    weight = np.where(np.abs(eta)<1.3, cenweight, forweight)
    return genw*weight



def get_ecal_bad_calib(run_number, lumi_number, event_number, year, dataset):
    bad = {}
    bad["2016"] = {}
    bad["2017"] = {}
    bad["2018"] = {}
    bad["2016"]["MET"]            = "ecalBadCalib/Run2016_MET.root"
    bad["2016"]["SinglePhoton"]   = "ecalBadCalib/Run2016_SinglePhoton.root"
    bad["2016"]["SingleElectron"] = "ecalBadCalib/Run2016_SingleElectron.root"
    bad["2017"]["MET"]            = "ecalBadCalib/Run2017_MET.root"
    bad["2017"]["SinglePhoton"]   = "ecalBadCalib/Run2017_SinglePhoton.root"
    bad["2017"]["SingleElectron"] = "ecalBadCalib/Run2017_SingleElectron.root"
    bad["2018"]["MET"]            = "ecalBadCalib/Run2018_MET.root"
    bad["2018"]["EGamma"]         = "ecalBadCalib/Run2018_EGamma.root"
    
    regular_dataset = ""
    regular_dataset = [name for name in ["MET","SinglePhoton","SingleElectron","EGamma"] if (name in dataset)]
    fbad = uproot3.open(bad[year][regular_dataset[0]])
    bad_tree = fbad["vetoEvents"]
    runs_to_veto = bad_tree.array("Run")
    lumis_to_veto = bad_tree.array("LS")
    events_to_veto = bad_tree.array("Event")

    # We want events that do NOT have (a vetoed run AND a vetoed LS and a vetoed event number)
    return np.logical_not(np.isin(run_number, runs_to_veto) * np.isin(lumi_number, lumis_to_veto) * np.isin(event_number, events_to_veto))


############
## BTag
## Btag recommend: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation#Recommendation_for_13_TeV_Data



from coffea.lookup_tools.correctionlib_wrapper import correctionlib_wrapper
from coffea.lookup_tools.dense_lookup import dense_lookup

class BTagCorrector:

    def __init__(self, tagger, year, workingpoint):
        self._year = year
        common = load('data/common.coffea')
        self._wp = common['btagWPs'][tagger][year][workingpoint]
        
        btvjson = correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'_UL/btagging.json.gz')
        self.sf = btvjson # ["deepJet_comb", "deepCSV_comb"]

        files = {
            '2016preVFP': 'btageff2016.merged',
            '2016postVFP': 'btageff2016.merged',
            '2017': 'btageff2017.merged',
            '2018': 'btageff2018.merged',
        }
        filename = 'hists/'+files[year]
        btag = load(filename)
        bpass = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag', 'pass').values()[()]
        ball = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag').values()[()]
        ball[ball<=0.]=1.
        nom = bpass / np.maximum(ball, 1.)
        self.eff = lookup_tools.dense_lookup.dense_lookup(nom, [ax.edges() for ax in btag[tagger].axes()[3:]])

    def btag_weight(self, pt, eta, flavor, istag):
        abseta = abs(eta)
        
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1b_Event_reweighting_using_scale
        def P(eff):
            weight = eff.ones_like()
            weight[istag] = eff[istag]
            weight[~istag] = (1 - eff[~istag])
            return weight.prod()

        '''
        Correction deepJet_comb has 5 inputs
        Input systematic (string): 
        Input working_point (string): L/M/T
        Input flavor (int): hadron flavor definition: 5=b, 4=c, 0=udsg
        Input abseta (real):
        Input pt (real):
        '''

        bc = flavor > 0
        light = ~bc
        
        eff = self.eff(flavor, pt, abseta)
        
        #sf_nom = self.sf.eval('central', flavor, abseta, pt)
        sf_nom = self.sf["deepJet_comb"].evaluate('central','M', flavor, abseta, pt)
        
        bc_sf_up_correlated = pt.ones_like()
        bc_sf_up_correlated[~bc] = sf_nom[~bc]
        bc_sf_up_correlated[bc] = self.sf["deepJet_comb"].evaluate('up_correlated', 'M', flavor, eta, pt)[bc]
        
        bc_sf_down_correlated = pt.ones_like()
        bc_sf_down_correlated[~bc] = sf_nom[~bc]
        bc_sf_down_correlated[bc] = self.sf["deepJet_comb"].evaluate('down_correlated', 'M', flavor, eta, pt)[bc]

        bc_sf_up_uncorrelated = pt.ones_like()
        bc_sf_up_uncorrelated[~bc] = sf_nom[~bc]
        bc_sf_up_uncorrelated[bc] = self.sf["deepJet_comb"].evaluate('up_uncorrelated', 'M', flavor, eta, pt)[bc]

        bc_sf_down_uncorrelated = pt.ones_like()
        bc_sf_down_uncorrelated[~bc] = sf_nom[~bc]
        bc_sf_down_uncorrelated[bc] = self.sf["deepJet_comb"].evaluate('down_uncorrelated', 'M', flavor, eta, pt)[bc]

        light_sf_up_correlated = pt.ones_like()
        light_sf_up_correlated[~light] = sf_nom[~light]
        light_sf_up_correlated[light] = self.sf["deepJet_comb"].evaluate('up_correlated', 'M', flavor, abseta, pt)[light]

        light_sf_down_correlated = pt.ones_like()
        light_sf_down_correlated[~light] = sf_nom[~light]
        light_sf_down_correlated[light] = self.sf["deepJet_comb"].evaluate('down_correlated', 'M', flavor, abseta, pt)[light]

        light_sf_up_uncorrelated = pt.ones_like()
        light_sf_up_uncorrelated[~light] = sf_nom[~light]
        light_sf_up_uncorrelated[light] = self.sf["deepJet_comb"].evaluate('up_uncorrelated', 'M', flavor, abseta, pt)[light]

        light_sf_down_uncorrelated = pt.ones_like()
        light_sf_down_uncorrelated[~light] = sf_nom[~light]
        light_sf_down_uncorrelated[light] = self.sf["deepJet_comb"].evaluate('down_uncorrelated', 'M', flavor, abseta, pt)[light]



        eff_data_nom  = np.minimum(1., sf_nom*eff)
        bc_eff_data_up_correlated   = np.minimum(1., bc_sf_up_correlated*eff)
        bc_eff_data_down_correlated = np.minimum(1., bc_sf_down_correlated*eff)
        bc_eff_data_up_uncorrelated   = np.minimum(1., bc_sf_up_uncorrelated*eff)
        bc_eff_data_down_uncorrelated = np.minimum(1., bc_sf_down_uncorrelated*eff)
        light_eff_data_up_correlated   = np.minimum(1., light_sf_up_correlated*eff)
        light_eff_data_down_correlated = np.minimum(1., light_sf_down_correlated*eff)
        light_eff_data_up_uncorrelated   = np.minimum(1., light_sf_up_uncorrelated*eff)
        light_eff_data_down_uncorrelated = np.minimum(1., light_sf_down_uncorrelated*eff)
       
        nom = P(eff_data_nom)/P(eff)
        bc_up_correlated = P(bc_eff_data_up_correlated)/P(eff)
        bc_down_correlated = P(bc_eff_data_down_correlated)/P(eff)
        bc_up_uncorrelated = P(bc_eff_data_up_uncorrelated)/P(eff)
        bc_down_uncorrelated = P(bc_eff_data_down_uncorrelated)/P(eff)
        light_up_correlated = P(light_eff_data_up_correlated)/P(eff)
        light_down_correlated = P(light_eff_data_down_correlated)/P(eff)
        light_up_uncorrelated = P(light_eff_data_up_uncorrelated)/P(eff)
        light_down_uncorrelated = P(light_eff_data_down_uncorrelated)/P(eff)
        '''
        print('nom',sf_nom)
        print('bc_up_correlated',bc_sf_up_correlated)
        print('bc_down_correlated',bc_sf_down_correlated)
        print('bc_up_uncorrelated',bc_sf_up_uncorrelated)
        print('bc_down_uncorrelated',bc_sf_down_uncorrelated)
        print('light_up_correlated',light_sf_up_correlated)
        print('light_down_correlated',light_sf_down_correlated)
        print('light_up_uncorrelated',light_sf_up_uncorrelated)
        print('light_down_uncorrelated',light_sf_down_uncorrelated)
        '''
        return np.nan_to_num(nom, nan=1.), np.nan_to_num(bc_up_correlated, nan=1.), np.nan_to_num(bc_down_correlated, nan=1.), np.nan_to_num(bc_up_uncorrelated, nan=1.), np.nan_to_num(bc_down_uncorrelated, nan=1.), np.nan_to_num(light_up_correlated, nan=1.), np.nan_to_num(light_down_correlated, nan=1.), np.nan_to_num(light_up_uncorrelated, nan=1.), np.nan_to_num(light_down_uncorrelated, nan=1.)


get_btag_weight = {
    'deepflav': {
        '2016preVFP': {
            'loose'  : BTagCorrector('deepflav','2016preVFP','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2016preVFP','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2016preVFP','tight').btag_weight
        },
        '2016postVFP': {
            'loose'  : BTagCorrector('deepflav','2016postVFP','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2016postVFP','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2016postVFP','tight').btag_weight
        },
        '2017': {
            'loose'  : BTagCorrector('deepflav','2017','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2017','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2017','tight').btag_weight
        },
        '2018': {
            'loose'  : BTagCorrector('deepflav','2018','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2018','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2018','tight').btag_weight
        }
    },
    'deepcsv' : {
        '2016preVFP': {
            'loose'  : BTagCorrector('deepcsv','2016preVFP','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2016preVFP','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2016preVFP','tight').btag_weight
        },
        '2016postVFP': {
            'loose'  : BTagCorrector('deepcsv','2016postVFP','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2016postVFP','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2016postVFP','tight').btag_weight
        },
        '2017': {
            'loose'  : BTagCorrector('deepcsv','2017','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2017','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2017','tight').btag_weight
        },
        '2018': {
            'loose'  : BTagCorrector('deepcsv','2018','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2018','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2018','tight').btag_weight
        }
    }
}


class DoubleBTagCorrector:

    def __init__(self, year):
        self._year = year
        sf = {
            '2018': {
                'value': np.array([0.82, 0.82, 0.75, 0.81]),
                'unc': np.array([np.sqrt(0.07**2 + 0.11**2), np.sqrt(0.07**2 + 0.11**2), np.sqrt(0.06**2 + 0.06**2), np.sqrt(0.05**2 + 0.01**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
            '2017': {
                'value': np.array([0.84, 0.84, 0.98, 0.86]),
                'unc': np.array([np.sqrt(0.05**2 + 0.13**2), np.sqrt(0.05**2 + 0.13**2), np.sqrt(0.05**2 + 0.12**2), np.sqrt(0.05**2 + 0.05**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
            '2016': {
                'value': np.array([1.01, 1.01, 0.95, 0.99]),
                'unc': np.array([np.sqrt(0.06**2 + 0.02**2), np.sqrt(0.06**2 + 0.02**2), np.sqrt(0.05**2 + 0.09**2), np.sqrt(0.06**2 + 0.00**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
        }
        self.sf_nom=lookup_tools.dense_lookup.dense_lookup(sf[year]['value'], sf[year]['edges'])
        self.sf_up=lookup_tools.dense_lookup.dense_lookup(sf[year]['value']+sf[year]['unc'], sf[year]['edges'])
        self.sf_down=lookup_tools.dense_lookup.dense_lookup(sf[year]['value']-sf[year]['unc'], sf[year]['edges'])
        

    def doublebtag_weight(self, pt):
        return  self.sf_nom(pt), self.sf_up(pt), self.sf_down(pt)

get_doublebtag_weight = {
    '2016': DoubleBTagCorrector('2016').doublebtag_weight,
    '2017': DoubleBTagCorrector('2017').doublebtag_weight,
    '2018': DoubleBTagCorrector('2018').doublebtag_weight,
}

class Reweighting:

    def __init__(self, year):
        self._year = year
        files = {
            '2016': 'reweighting2016.scaled',
            '2017': 'reweighting2017.scaled',
            '2018': 'reweighting2018.scaled',
        }
        filename = 'hists/'+files[year]
        reweighting = load(filename)
        hists = load(filename)
        if year == '2018':
            hists['data']['reweighting'].scale(2.)
        num=hists['data']['reweighting'].integrate('process').values()[()]
        den=hists['bkg']['reweighting'].integrate('process').values()[()]
        reweighting = num / np.maximum(den, 1.)
        self.reweighting = lookup_tools.dense_lookup.dense_lookup(reweighting, [ax.edges() for ax in hists['data']['reweighting'].axes()[1:]])

    def weight(self, tau21, pt, eta):
        weight = self.reweighting(tau21, pt, eta)
        return weight

get_reweighting = {
    '2016': Reweighting('2016').weight,
    '2017': Reweighting('2017').weight,
    '2018': Reweighting('2018').weight,
}


####
# Jet
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
# But it still uses txt file
####

Jetext = extractor()
for directory in ['jec_UL', 'jersf_UL', 'jr_UL', 'junc_UL']:
    directory='data/'+directory
    print('Loading files in:',directory)
    for filename in os.listdir(directory):
        if '~' in filename: continue
#        if 'DATA' in filename: continue
        if "Regrouped" in filename: continue
        if "UncertaintySources" in filename: continue
        if 'AK4PFchs' in filename:
            filename=directory+'/'+filename
            print('Loading file:',filename)
            Jetext.add_weight_sets(['* * '+filename])
        if 'AK8' in filename:
            filename=directory+'/'+filename
            print('Loading file:',filename)
            Jetext.add_weight_sets(['* * '+filename])
    print('All files in',directory,'loaded')
Jetext.finalize()
Jetevaluator = Jetext.make_evaluator()


corrections = {}
corrections = {
    'get_met_trig_weight':      get_met_trig_weight,
    'get_ele_loose_id_sf':      get_ele_loose_id_sf,
    'get_ele_tight_id_sf':      get_ele_tight_id_sf,
    'get_ele_trig_weight':      get_ele_trig_weight,
    'get_ele_reco_sf_below20':  get_ele_reco_sf_below20,
    'get_ele_reco_err_below20': get_ele_reco_err_below20,
    'get_ele_reco_sf_above20':  get_ele_reco_sf_above20,
    'get_ele_reco_err_above20': get_ele_reco_err_above20,
    'get_pho_loose_id_sf':      get_pho_loose_id_sf,
    'get_pho_tight_id_sf':      get_pho_tight_id_sf,
    'get_pho_trig_weight':      get_pho_trig_weight,
    'XY_MET_Correction':        XY_MET_Correction,
    'pu_weight':                pu_weight,
    'get_nlo_ewk_weight':       get_nlo_ewk_weight,
    'get_ttbar_weight':         get_ttbar_weight,
    'get_msd_weight':           get_msd_weight,
    'get_btag_weight':          get_btag_weight,
    'jetevaluator':             Jetevaluator,
}


save(corrections, 'data/testcorrections.coffea')

