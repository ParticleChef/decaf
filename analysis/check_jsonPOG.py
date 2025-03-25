import numpy as np
import os
import uproot
from coffea import hist, nanoevents, util
from coffea.util import load, save
import coffea.processor as processor
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import correctionlib
from coffea.nanoevents.methods import vector

def check_json(path):
    year = '2022pre'
    evaluator = correctionlib.CorrectionSet.from_file(path)
    for corr in evaluator.values():
        if not corr.name == "NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": continue
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")


check_json("/u/user/particlechef/nanoaod_study/decaf/analysis/data/MuonSF/2023pre/muon_Z.json.gz")
def get_mu_highpt_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_HighPtID_DEN_TrackerMuons"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

def get_mu_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(flateta, flatpt, "nominal")

    return ak.unflatten(weight, counts=counts)

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


isRealsample = True
sample = {
    '18WJet': "root://dcache-cms-xrootd.desy.de:1094//store/mc/Run3Summer22NanoAODv12/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/055aacfc-52a3-4274-bfc6-767d6f193e9c.root",
}
sample_name = '18WJet'

if isRealsample:
    events = NanoEventsFactory.from_root(
            sample[sample_name],
            entry_stop=100_000,
            schemaclass=NanoAODSchema,
            ).events()

mu = events.Muon
#puweight = get_pu_weight('2022post', events.Pileup.nTrueInt )
#print(puweight)
#mu['isloose'] = isLooseMuon(mu, year)
#mu_loose = mu[mu.isloose]
for y in ['2022pre', '2022post', '2023pre', '2023post']:
    loose_mu_sf = get_mu_loose_id_sf(y, abs(mu.eta), mu.pt)
    print(loose_mu_sf)
