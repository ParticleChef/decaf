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


year = '2018'
isRealsample = False
sample = {
    'kit': "./selectedEvents_fromKIT.root",
    'ntp': "/10T/scratch/jhong/selectedEvents_result_files/merged_result.root",
    "GJet": "/data/mc/UL2018NanoAODv9//G1Jet_LHEGpT-150To250_TuneCP5_13TeV-amcatnlo-pythia8/MC_2018_NanoAODv9_1-32.root",
    "DY" : "/data/mc/UL2018NanoAODv9/DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/MC_2018_NanoAODv9_57.root",
    "1" : "./MC_2018_NanoAODv9_20.root"
}
sample_name = 'kit'

if isRealsample:
    events = NanoEventsFactory.from_root(
            sample[sample_name],
            entry_stop=100_000,
            schemaclass=NanoAODSchema,
            ).events()

events = uproot.open(sample[sample_name]+":Events")
#print(events.keys())
#print(events['LooseMuon_Pt'])
nevt = 1000

###
event_number = events['Evt_ID'].array()
test_evtID = ak.Array([1564487928, 1564201472, 1564036282, 206041131, 180664680, 261568833, 56316591])
selected_mask = np.isin(event_number,test_evtID)
event_number_new = event_number[selected_mask]
#print(len(event_number_new), len(test_evtID))
#print(selected_mask)
LooseMuon_pt = events['LooseMuon_Pt'].array()[selected_mask]
METpt  = events['MET_T1_Pt_nom'].array()[selected_mask]
METphi = events['MET_T1_Phi_nom'].array()[selected_mask]
recoilpt  = events['Hadr_Recoil_MET_T1_Pt_nom'].array()[selected_mask]
recoilphi = events['Hadr_Recoil_MET_T1XY_Phi_nom'].array()[selected_mask]
muon_pt = events['LooseMuon_Pt'].array()[selected_mask]
ele_pt = events['LooseElectron_Pt'].array()[selected_mask]
pho_pt = events['LoosePhoton_Pt'].array()[selected_mask]
muon_phi = events['LooseMuon_Phi'].array()[selected_mask]
ele_phi = events['LooseElectron_Phi'].array()[selected_mask]
pho_phi = events['LoosePhoton_Phi'].array()[selected_mask]


for i in range(len(METpt)):
    print('muon_pt: ', muon_pt[i],' phi: ',muon_phi[i], ' ID: ', muon.looseId[i], ' isoID: ', muon.pfRelIso04_all[i])
    print('elec_pt: ', ele_pt[i], ' phi: ',ele_phi[i], ' ID: ', ele.cutBased[i])
    print()


print(' METpt    METphi    ele_pt    ele_phi   pho_pt   pho_phi   mu_pt    mu_phi    recoil  ')
for i in range(len(METpt)):
    print('%6.2f' % METpt[i],'   %6.2f' % METphi[i], '   %6s' %ele_pt[i], '   %8s' %ele_phi[i], ' %6s' % pho_pt[i], ' %8s' % pho_phi[i],' %5s' % muon_pt[i],' % 8s' % muon_phi[i], '  %6.2f' %recoilpt[i])
    if len(ele_pt[i]) > 0:
        xxyy = (METpt[i]*np.cos(METphi[i])+ele_pt[i]*np.cos(ele_phi[i]))*(METpt[i]*np.cos(METphi[i])+ele_pt[i]*np.cos(ele_phi[i]))+(METpt[i]*np.sin(METphi[i])+ele_pt[i]*np.sin(ele_phi[i]))*(METpt[i]*np.sin(METphi[i])+ele_pt[i]*np.sin(ele_phi[i]))
        print('  ==> MET+ele pt: ',np.sqrt(xxyy))
    if len(muon_pt[i]) > 0:
        xxyy = (METpt[i]*np.cos(METphi[i])+muon_pt[i]*np.cos(muon_phi[i]))*(METpt[i]*np.cos(METphi[i])+muon_pt[i]*np.cos(muon_phi[i]))+(METpt[i]*np.sin(METphi[i])+muon_pt[i]*np.sin(muon_phi[i]))*(METpt[i]*np.sin(METphi[i])+muon_pt[i]*np.sin(muon_phi[i]))
        print('  ==> MET+muon pt: ',np.sqrt(xxyy))
    if len(pho_pt[i]) > 0:
        xxyy = (METpt[i]*np.cos(METphi[i])+pho_pt[i]*np.cos(pho_phi[i]))*(METpt[i]*np.cos(METphi[i])+pho_pt[i]*np.cos(pho_phi[i]))+(METpt[i]*np.sin(METphi[i])+pho_pt[i]*np.sin(pho_phi[i]))*(METpt[i]*np.sin(METphi[i])+pho_pt[i]*np.sin(pho_phi[i]))
        print('  ==> MET+pho pt: ',np.sqrt(xxyy))
    print()

#print('MET: ', METpt, ', ')
#print('e:   ', ele_pt, ', ')
#print('mu:  ', muon_pt, ', ')
#print('pho: ', pho_pt, ', ')
#print('recoil: ', recoilpt, ', ', recoilphi)
###
test_evtID = [1564487928, 1564201472, 1564036282, 206041131, 180664680, 261568833, 56316591] 
for i in range(len(event_number)):
#for i in test_evtID:
    if event_number[i] in test_evtID:
        continue
        print('event number: ', event_number[i])
        noflep = len(events['LooseElectron_Pt'].array()[i])+len(events['LoosePhoton_Pt'].array()[i])+len(events['LooseMuon_Pt'].array()[i])
        if noflep > 0:
            print('   %6s' %events['LooseElectron_Pt'].array()[i], ' %6s' % events['LoosePhoton_Pt'].array()[i], ' %5s' % events['LooseMuon_Pt'].array()[i])
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
#e = events.Electron
#e['isloose'] = isLooseElectron(e, year)
#e['T'] = ak.zip({ "r": e.pt, "phi": e.phi,},
#    with_name="PolarTwoVector",
#    behavior=vector.behavior,
#)
#e_loose = e[e.isloose]
#
#pho = events.Photon
#pho['T'] = ak.zip({ "r": pho.pt, "phi": pho.phi,},
#    with_name="PolarTwoVector",
#    behavior=vector.behavior,
#)
#pho_tight = pho[isTightPho(pho,'2018')]
#pho_loose = pho[isLoosePho(pho,'2018')]

#for i in range(nevt):
#    #if len(mu_loose[i]) == 0: continue
#    if not events.event[i] in test_evtID: continue
#    print('event ', events.event[i])
#    print('met      : %.3f' % met.pt[i])
#    print('loose mu : ', mu_loose.pt[i])
#    print('loose e  : ', e_loose.pt[i])
#    print('loose pho: ', pho_loose.pt[i])
#    temp_vector = e_loose.T[i] + pho_loose.T[i]
#    met_lep = met[i] + mu_loose.T[i] + temp_vector #e_loose.T[i] + pho_loose.T[i]
#    print('met + lep: %s' %met_lep)
#    print()


def testRecoil():
    recoil = met + pho
    a = 33.7*np.sin(1.19)+191*np.sin(-2.34)
    b = 33.7*np.cos(1.19)+191*np.cos(-2.34)
    print("met pt, phi: ", met.pt, met.phi)
    print("pho pt, phi: ", pho.pt, pho.phi)
    print("rec pt, phi: ", recoil.pt, recoil.phi)
    print("a and b ", a, b)
    print("vec sum ", recoil)

ni = 0
nf = len(events)

def rochesterTest():
    get_rochester = corrections['get_mu_rochester_sf']['2018']
    
    mu = events.Muon
    kspread = get_rochester.kSpreadMC(mu.charge, mu.pt, mu.eta, mu.phi, mu.matched_gen.pt)
    mc_rand = ak.unflatten(np.random.rand(ak.count(ak.flatten(mu.pt))), ak.num(mu.pt))
    ksmear = get_rochester.kSmearMC(mu.charge, mu.pt, mu.eta, mu.phi, mu.nTrackerLayers, mc_rand)
    hasgen = ~np.isnan(ak.fill_none(events.Muon.matched_gen.pt, np.nan))
    k = ak.where(hasgen, kspread, ksmear)
    
    print('kspread ', kspread)
    print('ksmear ', ksmear)
    print('rochester k = ', k)
    isData = False
    if isData:
        k = get_mu_rochester_sf.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi)

def print_4_line(photon, ievt):
    print("%20s" % photon.pt[ievt], "%20s" % photon.phi[ievt], "%20s" % photon.eta[ievt], "%20s" % photon.cutBased[ievt], end='')

def check_pho():
    title_flag = True
    for i in range(ni, nf):
        if title_flag:
            print("%60s" % str("Original"), "%60s" % str("Tight"), "%60s" % str("Loose"))
            for i in range(3):
                print("%20s" % str("pt"), "%20s" % str("phi"), "%20s" % str("eta"),  "%20s" % str("cutBased"), end='')
            print()
            title_flag = False
        if len(pho_tight[i]) > 0 and len(pho_loose[i]) == 0:
            print_4_line(pho, i)
            print_4_line(pho_tight, i)
            print_4_line(pho_loose,i)
            print()

def print_pho(photon, name, nevti, nevtf):
    print(str(name))
    print('pt         : ', end=' ')
    for i in range(nevti, nevtf):
        print(photon.pt[i], end=' ')
    print()
    print('phi        : ', end=' ')
    for i in range(nevti, nevtf):
        print(photon.phi[i], end=' ')
    print()
    print('cutBased ID: ', end=' ')
    for i in range(nevti, nevtf):
        print(photon.cutBased[i], end=' ')
    print()
    print()
#print_pho(pho, 'Original Photon', ni,nf)
#print_pho(pho_tight, 'Tight Photon', ni,nf)
#print_pho(pho_loose, 'Loose Photon', ni,nf)

def def_objects():
    fj = events.AK15PFPuppi
    fj["pt"] = events.AK15PFPuppi.Jet_pt
    fj["eta"]  = events.AK15PFPuppi.Jet_eta
    fj["phi"] = events.AK15PFPuppi.Jet_phi
    fj['isgood'] = (fj.Jet_pt > 160) & (abs(fj.Jet_eta) < 2.4)
    fj_good = fj[fj.isgood]
    leading_fj = ak.firsts(fj_good)
    
    
    j = events.Jet
    j['isgood'] = (j.pt > 160) & (abs(j.eta) < 2.4)
    #j['isiso'] = ak.all(j.metric_table(leading_fj) > 1.5, axis=2)
    j_good = j[j.isgood]
    j_iso_mask = ak.all(j_good.metric_table(leading_fj) > 1.5, axis=2)
    dphi_j_fj = j_good.metric_table(leading_fj)
    dphi_j_fj2= j_good.delta_r(leading_fj)

def j_fj_test():
    print("fj pt: ", fj_good.Jet_pt)
    print("fj et: ", fj_good.Jet_eta)
    print()
    print(">>>>>>> PRINT FJ INFO <<<<<<<<")
    print("leading fj pt: ", leading_fj.Jet_pt)
    print("leading fj et: ", leading_fj.Jet_eta)
    print()
    print(">>>>>>> PRINT J INFO <<<<<<<<<")
    print("j pt: ", j_good.pt)
    print("j et: ", j_good.eta)
    print()
    print("j,leading_fj metric_table")
    print(" >>>>>> ", dphi_j_fj)
    print()
    print("j,leading_fj delta_phi")
    print(" >>>>>> ", dphi_j_fj2)
    print()
    print()

def ak_any_all():
    sample_array = np.array([[1,2],[0,0],[34,5],[50,55]])
    ak_all_array = ak.all(sample_array<6, axis=1)
    ak_any_array = ak.any(sample_array<6, axis=1)
    print(sample_array)
    print("Q: array < 6")
    print("ak.all ",ak_all_array)
    print("ak_any ",ak_any_array)



def check_mask(event, mask, loop, skip=False):
    event['ismask'] = mask
    masked_event = event[event.ismask]
    for i in range(loop):
        #if skip and ak.all(pho.pt[i] < 230, axis=0): continue
        #if skip and ak.all(event.ismask[i], axis=0): continue
        if skip and ak.all(mask==False, axis=1)[i]: continue
        print('%3d' % i, ' pt: ',masked_event.pt[i], ' eta: ', masked_event.eta[i], ' id: ', masked_event.cutBased[i])


def id_test():
    print()
    print('=== all ')
    
    
    print()
    print('=== pt > 230 ')
    check_mask(pho, pho.pt > 230, nevt, True)
    print()
    print('=== tight id')
    check_mask(pho, pho.cutBased==2, nevt, True)
    
    mask = ~np.isnan(ak.ones_like(pt))
    mask = (pt > 230) & (tight_id == 2) & ( abs(eta) < 1.479)
    tiphomask = mask&(pho.electronVeto)&(pho.isScEtaEB) #tight photons are barrel only
    
    
    print()
    print('=== eta < 1.479')
    check_mask(pho, abs(pho.eta) < 1.479, nevt, True)
    
    print()
    print('========== tight photon id ==========')
    print('== pt > 230, |eta| < 1.479, tight-id == 2')
    check_mask(pho, mask, nevt, True)



print()


def and_test():
    sample = np.array([[1,2,3],[0,0],[34,5],[50,55,1]])
    a = np.array([True,True,False, True, True,False])
    b = np.array([True,False,False,True,False,False])
    c = np.array([True,True,True,False,False,False])
    
    print('a         ', a)
    print('b         ', b)
    print('c         ', c)
    print('a & b     ', a&b)
    print('a & b & c ', a&b&c)




def check_json():
    year = '2018'
    evaluator = correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'_UL/btagging.json.gz')
    for corr in evaluator.values():
        if not corr.name == "deepJet_comb": continue
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")


#j, nj = ak.flatten(j), ak.num(j)

#j_flv = np.array(j.hadronFlavour)
#j_eta = np.array(abs(j.eta))
#j_pt  = np.array(j.pt)



