#!/usr/bin/env python
import numpy as np
import awkward as ak
import json
from correctionlib.schemav2 import CorrectionSet
import gzip
import correctionlib
from correctionlib import convert

aaa = '''
year = ["2022pre", "2022post", "2023pre", "2023post"]
cutbased = ["Loose", "Tight"]

for yr in year:
    for i in cutbased:
        trghist = 'data/ElectronTrigEff/'+yr+'/'+i+'/egammaEffi.txt_EGM2D.root:EGamma_SF2D'

        corr = convert.from_uproot_THx(trghist)
        cset = CorrectionSet(schema_version=2, corrections=[corr])

        with gzip.open("data/ElectronTrigEff/"+yr+"/"+i+"/egammaEffi"+i+yr+".json.gz", "wt") as f:
            f.write(cset.json(exclude_unset=True))
'''
## One loop
#trghist = 'data/ElectronTrigEff/'+year+'/'+cutbased+'/egammaEffi.txt_EGM2D.root:EGamma_SF2D'
#
#corr = convert.from_uproot_THx(trghist)
#cset = CorrectionSet(schema_version=2, corrections=[corr])
#
#with gzip.open("data/ElectronTrigEff/"+year+"/"+cutbased+"/egammaEffi"+cutbased+year+".json.gz", "wt") as f:
#    f.write(cset.json(exclude_unset=True))

pt  = np.array([25.1, 45, 601, 50, 51])
eta = np.array([1.2, 1.1, -2.1, 3.1, -3.1])


year = "2022pre"
cutbased = "Tight"

evaluator = correctionlib.CorrectionSet.from_file("data/ElectronTrigEff/"+year+"/"+cutbased+"/egammaEffi"+cutbased+year+".json.gz")

pt = ak.where((pt<25.0), ak.full_like(pt,25.), pt)
pt = ak.where((pt>500.0), ak.full_like(pt,499.9), pt)
print(pt)
eta = ak.where((eta>2.5), ak.full_like(eta,2.499), eta)
eta = ak.where((eta<-2.5), ak.full_like(eta,-2.499), eta)
print(eta)

#flatpt = ak.flatten(pt)
#flateta, counts = ak.flatten(eta), ak.num(eta)

sf = evaluator["h2_scaleFactorsEGamma"].evaluate(eta, pt)
#sf = evaluator["h2_scaleFactorsEGamma"].evaluate(flateta, flatpt)

print(sf)

