import numpy as np
import awkward as ak
import correctionlib
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

#Def get_nnlo_nlo_weight(correctionname, pt):
#    evaluator = correctionlib.CorrectionSet.from_file('data/VJets_SFs/kfactors_wjets.json.gz')
#    #evaluator = CorrectionSet.from_file('data/VJets_SFs/kfactors_wjets.json')
#
#    print(evaluator[correctionname].evaluate(pt))
#
#
#Get_nnlo_nlo_weight('electron_kfactor', 850)

year = '2022'
isRealsample = True
sample = {
    "1" : "/data/mc/privatemc/2022/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3_AK15_ParTv2_Run3Summer22MiniAODv4-130X_v5_ext1-v2/0000/nano_1.root"
}
sample_name = '1'

if isRealsample:
    events = NanoEventsFactory.from_root(
            sample[sample_name],
            entry_stop=100_000,
            schemaclass=NanoAODSchema,
            ).events()

print(len(events))
def get_nnlo_nlo_wjet(channel, mass):
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

    #for low, high, value in kfactors[channel]:
    #    if low <= mass < high:
    #        return value
    #raise ValueError(f"Mass {mass} GeV out of range for {channel} channel.")



gen = events.GenPart
gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
genWs = gen[gen.isW]

condition = True
print('first genWs', genWs.pt)
#nnlo_nlo = ak.where(
#    ((ak.num(genWs, axis=1)>0)&(ak.firsts(genWs).pt>=100)),
#    get_nnlo_nlo_wjet('electron', ak.firsts(genWs).pt),
#    np.ones(len(events), dtype='float')
#    )
nnlo_nlo = ak.where(
    ((ak.num(genWs, axis=1)>0)&(genWs.pt>=100)),
    get_nnlo_nlo_wjet('electron', genWs.pt),
    np.ones(len(events), dtype='float')
    )
print('sf', nnlo_nlo)
sf_list = {}
for sf in ak.flatten(nnlo_nlo):
    if not str(sf) in sf_list.keys():
        sf_list[str(sf)] = 1
    else:
        sf_list[str(sf)] = sf_list[str(sf)] + 1

print(sf_list)


# /data/mc/privatemc/2022/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3_AK15_ParTv2_Run3Summer22MiniAODv4-130X_v5_ext1-v2/0000/  nano_1.root
