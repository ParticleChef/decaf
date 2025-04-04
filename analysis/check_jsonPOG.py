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
        #if not corr.name == "NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": continue
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")


check_json("/home/jhong/nanoaod-study/decaf/analysis/data/MuonSF/2023pre/muon_Z.json.gz")

def get_electron_id_sf (year, wp, eta, pt, phi):
    evaluator = correctionlib.CorrectionSet.from_file('data/EGammaSF/'+year+'/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatpt = ak.flatten(pt)

    flatphi = ak.flatten(phi)
    yr = {
        '2022pre' : '2022Re-recoBCD',
        '2022post': '2022Re-recoE+PromptFG',
        '2023pre' : '2023PromptC',
        '2023post': '2023PromptD'
    }
    print(yr[year])
    print(wp)
    print(flateta)
    print(flatpt)
    print(flatphi)
    
    if '2022' in year:
        weight = evaluator["Electron-ID-SF"].evaluate(yr, "sf", wp, flateta, flatpt)
    if '2023' in year:
        weight = evaluator["Electron-ID-SF"].evaluate(yr, "sf", wp, flateta, flatpt, flatphi)

    return ak.unflatten(weight, counts=counts)

isRealsample = True
sample = {
    'TT_run3': "/data/mc/Run3Summer22NanoAODv12/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/f7267ea1-1288-4112-a06a-849bf1f36dfa.root"
#    '18WJet': "root://dcache-cms-xrootd.desy.de:1094//store/mc/Run3Summer22NanoAODv12/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/055aacfc-52a3-4274-bfc6-767d6f193e9c.root",
}
sample_name = 'TT_run3'

if isRealsample:
    events = NanoEventsFactory.from_root(
            sample[sample_name],
            entry_stop=100_000,
            schemaclass=NanoAODSchema,
            ).events()

ele = events.Photon
#puweight = get_pu_weight('2022post', events.Pileup.nTrueInt )
#print(puweight)
#mu['isloose'] = isLooseMuon(mu, year)
#mu_loose = mu[mu.isloose]
for y in ['2022pre', '2022post', '2023pre', '2023post']:
    loose_electron_sf = get_electron_id_sf(y, 'Loose', abs(ele.eta), ele.pt, ele.phi)
    print(loose_electron_sf)
