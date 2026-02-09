import numpy as np
import awkward as ak
from coffea.util import save
from coffea.nanoevents.methods import vector as v
from coffea import lookup_tools, jetmet_tools, util
from coffea.lookup_tools import extractor, dense_lookup
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory
import correctionlib
from correctionlib import convert

def isGoodJet(jet, year):
    pt = jet.pt
    eta = jet.eta
    chHEF = jet.chHEF
    neHEF = jet.neHEF
    chEmEF = jet.chEmEF
    neEmEF = jet.neEmEF
    muEF = jet.muEF
    chMultiplicity = jet.chMultiplicity
    neMultiplicity = jet.neMultiplicity
    multiplicity = chMultiplicity + neMultiplicity
    def getJetID(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity):
        evaluator = correctionlib.CorrectionSet.from_file('data/JMESF/'+year+'/jetid.json.gz')
        corr = evaluator["AK4PUPPI_TightLeptonVeto"]
        args = (
            float(eta),
            float(chHEF), float(neHEF), float(chEmEF), float(neEmEF), float(muEF),
            int(chMultiplicity), int(neMultiplicity), int(multiplicity),
        )
        out = corr.evaluate(*args)
        return out
    jetId = getJetID(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
    print(jetId)
    mask = (pt > 30) & (abs(eta) < 2.4) & (jetId == 1) & (jet.puId >= 4)
    return mask

pt = ak.Array([35, 40, 50])
eta = ak.Array([0.5, -1.2, 2.3])
chHEF = ak.Array([0.4, 0.3, 0.2])
neHEF = ak.Array([0.25, 0.35, 0.15])
chEmEF = ak.Array([0.1, 0.2, 0.15])
neEmEF = ak.Array([0.15, 0.1, 0.25])
muEF = ak.Array([0.1, 0.05, 0.2])
chMultiplicity = ak.Array([1, 2, 3])
neMultiplicity = ak.Array([2, 1, 0])
multiplicity = chMultiplicity + neMultiplicity
puId = ak.Array([7, 5, 6])  # 예시 puId 값
year = '2024'

jet = ak.zip({
    'pt': pt,
    'eta': eta,
    'chHEF': chHEF,
    'neHEF': neHEF,
    'chEmEF': chEmEF,
    'neEmEF': neEmEF,
    'muEF': muEF,
    'chMultiplicity': chMultiplicity,
    'neMultiplicity': neMultiplicity,
    'puId': puId
}, with_name="Momentum4D")
mask = isGoodJet(jet, year)
print("Good Jet Mask:", mask)