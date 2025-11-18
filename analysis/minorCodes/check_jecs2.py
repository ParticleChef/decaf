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






def jec_files(year, jec_algo, dojer):

    jec_level = ["L1FastJet", "L2L3Residual", "L2Relative", "L3Absolute", "Uncertainty"]
    jec_tag = {
        "2022pre"  : "Summer22_22Sep2023_V2_MC",
        "2022post" : "Summer22EE_22Sep2023_V2_MC",
        "2023pre"  : "Summer23Prompt23_V1_MC",
        "2023post" : "Summer23BPixPrompt23_V1_MC"
    }

    jer_level = ["PtResolution"]#,"SF"]
    jer_tag = {
        "2022pre"  : ["Summer22_22Sep2023_JRV1_MC"],
        "2022post" : ["Summer22EE_22Sep2023_JRV1_MC"],
        "2023pre"  : ["Summer23Prompt23_RunCv123_JRV1_MC", "Summer23Prompt23_RunCv4_JRV1_MC"],
        "2023post" : ["Summer23BPixPrompt23_RunD_JRV1_MC"]
    }

    filenames = []
    for lv in jec_level:
        filename = f"{jec_tag[year]}_{lv}_{jec_algo}.txt"
        if 'Unc' in lv: 
            filename = f"{jec_tag[year]}_{lv}_{jec_algo}.txt"
            #tempname = f"{jec_tag[year]}_{lv}_{jec_algo}.txt".split('_')
            #filename = tempname[0] + tempname[1] + '_' + '_' .join(tempname[2:])
            #print('filename unc: ', filename)
        filenames.append(filename)
    if dojer:
        for v in jer_tag[year]:
            filenames += [f"{v}_{lv}_{jec_algo}.txt" for lv in jer_level]

    return filenames

print(jec_files("2022pre", "AK4PFPuppi", False))


def jet_factory_factory(files, year):
    ext = extractor()
    directory='data/JetMETCorr/'+year
    for filename in files:
        ext.add_weight_sets([f"* * {directory+'/'+filename}"])
    ext.finalize()
    jec_stack = JECStack(ext.make_evaluator())
    return CorrectedJetsFactory(jec_name_map, jec_stack)

def jet_factory(year, jec_algo, dojer):
    return jet_factory_factory(jec_files(year, jec_algo, dojer), year)

def add_jec_variables(jets, event_rho):
	jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
	jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
	jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
	jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
	return jets



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



print("function jet_factory: ", jet_factory("2022pre","AK4PFPuppi", False))
jjjjj = jet_factory("2022pre","AK4PFPuppi", False)
print(events.Jet.pt)
jets = jjjjj.build(add_jec_variables(events.Jet, events.Rho.fixedGridRhoFastjetAll), jec_cache)
print(jets.pt)



