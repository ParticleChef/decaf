import numpy as np
import os
import uproot
import cachetools
from coffea import hist, nanoevents, util
from coffea.util import load, save
import coffea.processor as processor
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import correctionlib
from coffea.nanoevents.methods import vector
from coffea import lookup_tools, jetmet_tools, util
from coffea.lookup_tools import extractor, dense_lookup
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory

def check_json(path):
    year = '2022pre'
    evaluator = correctionlib.CorrectionSet.from_file(path)
    for corr in evaluator.values():
        #if not corr.name == "NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": continue
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        #print(f"{corr.name}")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")


check_json("./data/JetMETCorr/2022pre/met_xyCorrections_2022_2022.json.gz")

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






isRealsample = True
sample = {
    'TT_run3': "/data/mc/Run3Summer22NanoAODv12/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/f7267ea1-1288-4112-a06a-849bf1f36dfa.root"
#    '18WJet': "root://dcache-cms-xrootd.desy.de:1094//store/mc/Run3Summer22NanoAODv12/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/055aacfc-52a3-4274-bfc6-767d6f193e9c.root",
}
sample_name = 'TT_run3'
jec_cache = cachetools.Cache(np.inf)
if isRealsample:
    events = NanoEventsFactory.from_root(
            sample[sample_name],
            entry_stop=100_000,
            schemaclass=NanoAODSchema,
            ).events()

    ele = events.Photon
    #print("Before: ", events.Jet['pt'])
    #jets = jet_factory['2022premc'].build(add_jec_variables(events.Jet, events.Rho.fixedGridRhoFastjetAll), jec_cache)
    #jets = jet_factory('2022pre', 'AK4PFPuppi', True).build(add_jec_variables(events.Jet, events.Rho.fixedGridRhoFastjetAll), jec_cache)
    #print("After: ", jets['pt'])
    #ele_reco = get_ele_reco_sf_20to75('2022pre', ele.eta, ele.pt, ele.phi)
    #puweight = get_pu_weight('2022post', events.Pileup.nTrueInt )
    #print(puweight)
    #mu['isloose'] = isLooseMuon(mu, year)
    #mu_loose = mu[mu.isloose]
    for y in ['2022pre', '2022post', '2023pre', '2023post']:
        #test_sf = get_electron_id_sf(y,'Loose', abs(ele.eta), ele.pt, ele.phi)
        #test_sf = jec_AK4(y, abs(ele.eta), ele.pt)
        test_sf = get_ele_reco_sf_Above75(y, ele.eta, ele.pt, ele.phi)
        print(y)
        print(test_sf)
        #break

def get_ele_loose_id_sf (year, eta, pt, phi):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatphi, flatpt = ak.flatten(phi), ak.flatten(pt)

    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD'
    }
    if '2022' in year:
        weight = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Loose", flateta, flatpt)
    if '2023' in year:
        weight = evaluator["Electron-ID-SF"].evaluate(yr[year], "sf", "Loose", flateta, flatpt, flatphi)

    return ak.unflatten(weight, counts=counts)


#def jec_AK4(year, eta, pt):
#    evaluator = correctionlib.CorrectionSet.from_file("data/JetMETCorr/"+year+"/jet_jerc.json.gz")
#
#    flateta, counts = ak.flatten(eta), ak.num(eta)
#    flatpt = ak.flatten(pt)
#
#    jec_level = ["L1FastJet", "L2L3Residual", "L2L3Residual", "L3Absolute"]
#    jec_tag = {
#        "2022pre"  : "Summer22_22Sep2023_V2_MC",
#        "2022post" : "Summer22EE_22Sep2023_V2_MC",
#        "2023pre"  : "Summer23Prompt23_V1_MC",
#        "2023post" : "Summer23BPixPrompt23_V1_MC"
#    }
#
#    sf = evaluator[f"{jec_tag[year]}_{jec_level[0]}_AK4PFPuppi"].evaluate(eta, pt)
#
#    return ak.unflatten(sf, counts=counts)


