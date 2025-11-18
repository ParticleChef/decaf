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


#year = '2018'
#isRealsample = True
#sample = {
#    'ntp': "/10T/scratch/jhong/selectedEvents_result_files/merged_result.root",
#    "GJet": "/data/mc/UL2018NanoAODv9//G1Jet_LHEGpT-150To250_TuneCP5_13TeV-amcatnlo-pythia8/MC_2018_NanoAODv9_1-32.root",
#    "DY" : "/data/mc/UL2018NanoAODv9/DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/MC_2018_NanoAODv9_57.root",
#    "1" : "./MC_2018_NanoAODv9_20.root"
#}
#sample_name = 'ntp'
#
#if isRealsample:
#    events = NanoEventsFactory.from_root(
#            sample[sample_name],
#            entry_stop=100_000,
#            schemaclass=NanoAODSchema,
#            ).events()
#

def check_json():
    year = '2022pre'
    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'/muon_Z.json.gz')
    for corr in evaluator.values():
        #if not corr.name == "deepJet_comb": continue
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")


check_json()
#def get_mu_loose_id_sf (year, eta, pt):
#    evaluator = correctionlib.CorrectionSet.from_file('data/MuonSF/'+year+'_UL/muon_Z.json.gz')
#
#    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
#    flateta, counts = ak.flatten(eta), ak.num(eta)
#
#    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
#    flatpt = ak.flatten(pt)
#    
#    if year == '2022pre':
#        weight = evaluator["NUM_HighPtID_DEN_TrackerMuons"].evaluate(year+'_UL', flateta, flatpt, "sf")
#    else:
#        weight = evaluator["NUM_HighPtID_DEN_TrackerMuons"].evaluate(year+'_UL', flateta, flatpt, "sf")
#
#    return ak.unflatten(weight, counts=counts)
#ids = load('data/ids.coffea')
#isLooseMuon = ids['isLooseMuon']
#
#corrections = load('data/corrections.coffea')
#
#met = events.MET
#met['T'] = ak.zip({ "r": met.pt, "phi": met.phi,},
#    with_name="PolarTwoVector",
#    behavior=vector.behavior,
#)
#
#mu = events.Muon
#mu['isloose'] = isLooseMuon(mu, year)
#mu['T'] = ak.zip({ "r": mu.pt, "phi": mu.phi,},
#    with_name="PolarTwoVector",
#    behavior=vector.behavior,
#)
#mu_loose = mu[mu.isloose]
#
#
#
#def rochesterTest():
#    get_rochester = corrections['get_mu_rochester_sf']['2018']
#    
#    mu = events.Muon
#    kspread = get_rochester.kSpreadMC(mu.charge, mu.pt, mu.eta, mu.phi, mu.matched_gen.pt)
#    mc_rand = ak.unflatten(np.random.rand(ak.count(ak.flatten(mu.pt))), ak.num(mu.pt))
#    ksmear = get_rochester.kSmearMC(mu.charge, mu.pt, mu.eta, mu.phi, mu.nTrackerLayers, mc_rand)
#    hasgen = ~np.isnan(ak.fill_none(events.Muon.matched_gen.pt, np.nan))
#    k = ak.where(hasgen, kspread, ksmear)
#    
#    print('kspread ', kspread)
#    print('ksmear ', ksmear)
#    print('rochester k = ', k)
#    isData = False
#    if isData:
#        k = get_mu_rochester_sf.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi)
#
#
#def ak_any_all():
#    sample_array = np.array([[1,2],[0,0],[34,5],[50,55]])
#    ak_all_array = ak.all(sample_array<6, axis=1)
#    ak_any_array = ak.any(sample_array<6, axis=1)
#    print(sample_array)
#    print("Q: array < 6")
#    print("ak.all ",ak_all_array)
#    print("ak_any ",ak_any_array)
#
#
#
#def check_mask(event, mask, loop, skip=False):
#    event['ismask'] = mask
#    masked_event = event[event.ismask]
#    for i in range(loop):
#        #if skip and ak.all(pho.pt[i] < 230, axis=0): continue
#        #if skip and ak.all(event.ismask[i], axis=0): continue
#        if skip and ak.all(mask==False, axis=1)[i]: continue
#        print('%3d' % i, ' pt: ',masked_event.pt[i], ' eta: ', masked_event.eta[i], ' id: ', masked_event.cutBased[i])
#
#
#
#
#
#def and_test():
#    sample = np.array([[1,2,3],[0,0],[34,5],[50,55,1]])
#    a = np.array([True,True,False, True, True,False])
#    b = np.array([True,False,False,True,False,False])
#    c = np.array([True,True,True,False,False,False])
#    
#    print('a         ', a)
#    print('b         ', b)
#    print('c         ', c)
#    print('a & b     ', a&b)
#    print('a & b & c ', a&b&c)
#



#j, nj = ak.flatten(j), ak.num(j)

#j_flv = np.array(j.hadronFlavour)
#j_eta = np.array(abs(j.eta))
#j_pt  = np.array(j.pt)



