import numpy as np
import os
import uproot
from coffea import  nanoevents, util
from coffea.util import load, save
import coffea.processor as processor
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import correctionlib
from coffea.nanoevents.methods import vector
from correctionlib import CorrectionSet
from coffea.lookup_tools import extractor, dense_lookup
from coffea import lookup_tools, jetmet_tools, util
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory


year = '2022pre'
isRealsample = True
sample = {
    '2023': "/home/jhong/run3Monotop/decaf/analysis/temp_nano/nano_1.root",
    '2022': "/data/mc/privatemc/2022/WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3_AK15_ParTv2_Run3Summer22MiniAODv4-130X_v5-v3/0000/nano_11.root",
    '2022dat': "/data/data/privatedata/2022/EGamma/NanoTuples-AK15_ParTv2_Run2022D-22Sep2023-v1/251105_073349/0000/nano_823.root",
    'ntp': "/10T/scratch/jhong/selectedEvents_result_files/merged_result.root",
    "GJet": "/data/mc/UL2018NanoAODv9//G1Jet_LHEGpT-150To250_TuneCP5_13TeV-amcatnlo-pythia8/MC_2018_NanoAODv9_1-32.root",
    "DY" : "/data/mc/UL2018NanoAODv9/DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/MC_2018_NanoAODv9_57.root",
    "1" : "./MC_2018_NanoAODv9_20.root"
}
sample_name = '2022'

if isRealsample:
    events = NanoEventsFactory.from_root(
            sample[sample_name],
            entry_stop=100_000,
            schemaclass=NanoAODSchema,
            ).events()

def check_json(year):
    #year = '2022pre'
    evaluator = correctionlib.CorrectionSet.from_file('data/JetMETCorr/'+year+'/jet_jerc.json.gz')
    for corr in evaluator.values():
        if not "_L" in corr.name: continue
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")

#check_json()

corrections = load('data/corrections.coffea')
ids = load('data/ids.coffea')
isLooseMuon  = ids['isLooseMuon']
isTightMuon  = ids['isTightMuon']
get_met_trig_weight	  = corrections['get_met_trig_weight']
get_ele_trig_weight	  = corrections['get_ele_trig_weight']
get_pho_trig_weight	  = corrections['get_pho_trig_weight']
met = events.MET
mu = events.Muon
mu['isloose'] = isLooseMuon(mu,year)
mu['istight'] = isTightMuon(mu,year)
mu['T'] = ak.zip(
	{
		"r": mu.pt,
		"phi": mu.phi,
	},
	with_name="PolarTwoVector",
	behavior=vector.behavior,
)
mu_loose=mu[mu.isloose]
mu_tight=mu[mu.istight]

# define leading mu
leading_mu = ak.firsts(mu_tight)
u = met+leading_mu.T
met_sf, met_sf_up, met_sf_down = get_met_trig_weight(year, met.pt)
u_sf, u_sf_up, u_sf_down = get_met_trig_weight(year, u.r)
mask = (u.r > 250)
breaknumber = 0
number = 5
print('met pt  ', met.pt[number])
print('met phi ', met.phi[number])
print('mu pt ', leading_mu.pt[number])
print('mu phi', leading_mu.phi[number])
print('u ', u.r[number])
print('length', len(u))
for i in range(len(u)):
    if u.r[i] == None: continue
    if u.r[i] < 350: continue

    print('met', met.pt[i])
    print('ut ', u.pt[i])
    print('met sf', met_sf[i])
    print('u sf', u_sf[i])

    breaknumber = breaknumber +1

    #if breaknumber > 200: break

#def isGoodAK4(j, year):
#    
#    pt    = j.pt
#    eta   = j.eta
#    jet_id= j.jetId
#    #pu_id=j.puId
#    chHEF  = j.chHEF
#    neHEF  = j.neHEF
#    chEmEF = j.chEmEF
#    neEmEF = j.neEmEF
#    muEF   = j.muEF
#    chMultiplicity = j.chMultiplicity
#    neMultiplicity = j.neMultiplicity
#    multiplicity   = neMultiplicity + chMultiplicity
#
#    def getJetID(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity):
#        evaluator = correctionlib.CorrectionSet.from_file('data/JetMETCorr/'+year+'/jetid.json.gz')
#        corr = evaluator["AK4PUPPI_TightLeptonVeto"]
#        counts = ak.num(eta)
#        eta, chHEF, neHEF, chEmEF, neEmEF, muEF = ak.flatten(eta), ak.flatten(chHEF), ak.flatten(neHEF), ak.flatten(chEmEF), ak.flatten(neEmEF), ak.flatten(muEF)
#        chMultiplicity, neMultiplicity, multiplicity = ak.flatten(chMultiplicity), ak.flatten(neMultiplicity), ak.flatten(multiplicity)
#
#        print('1', eta,ak.type(eta))
#        print('2', chHEF,ak.type(chHEF))
#        print('3', neHEF,ak.type(neHEF))
#        print('4', chEmEF,ak.type(chEmEF))
#        print('5', neEmEF,ak.type(neEmEF))
#        print('6', muEF,ak.type(muEF))
#        print('7', chMultiplicity,ak.type(chMultiplicity))
#        print('8', neMultiplicity,ak.type(neMultiplicity))
#        print('9', multiplicity,ak.type(multiplicity))
#        args = (
#            eta,
#            chHEF, neHEF, chEmEF, neEmEF, muEF,
#            chMultiplicity, neMultiplicity, multiplicity
#        )
#        #out = corr.evaluate(*args)
#        out = corr.evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
#        return ak.unflatten(out, counts)
#    
#    jetId = getJetID(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
#    print('jetId', jetId)
#    mask = (pt > 30) & (abs(eta) < 2.4) & (jetId == 1)
#
#    return mask
#j = events.Jet
#isgood = isGoodAK4(j, '2022pre')
#print('isgood jet: ', isgood)

#get_ele_trig_weight      = corrections['get_ele_trig_weight']
#e = events.Electron
#e_weight = get_ele_trig_weight('2022pre', e.eta, e.pt, 'Loose')
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



