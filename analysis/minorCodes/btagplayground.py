#! /usr/bin/env python
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

from correctionlib import convert

import json

import hist

from coffea.lookup_tools.correctionlib_wrapper import correctionlib_wrapper
from coffea.lookup_tools.dense_lookup import dense_lookup

year = '2022pre' ## '2022post', '2023pre', '2023post'
isRealsample = True
sample = {
    '2022': "/data/mc/privatemc/2022/WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3_AK15_ParTv2_Run3Summer22MiniAODv4-130X_v5-v3/0000/nano_1.root"
}
sample_name = '2022'

if isRealsample:
    events = NanoEventsFactory.from_root(
            sample[sample_name],
            entry_stop=100_000,
            schemaclass=NanoAODSchema,
            ).events()

class BTagCorrector:

    def __init__(self, tagger, year, workingpoint):
        self._year = year

        wp = {}
        wp['loose'] = 'L'
        wp['medium'] = 'M'
        wp['tight'] = 'T'
        self._wp = wp[workingpoint]

        btvjson = {}
        if year == 2024:
            btvjson['PNetUParT'] = {
                'comb': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["UParTAK4_comb"],
                'ligh': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["UParTAK4_light"],
            }
        else:
            btvjson['PNetUParT'] = {
                'comb': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["particleNet_comb"],
                'ligh': correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'/btagging.json.gz')["particleNet_light"],
            }
        self.sf = btvjson[tagger]

        files = {
            '2022pre': 'btageff2022.merged',
        }
        filename = 'hists/'+files[year]
        btag_file = load(filename)
        for k in btag_file[tagger]:
            try:
                btag += btag_file[tagger][k]
            except:
                btag = btag_file[tagger][k]
        print('btag ', btag)

        bpass = btag[{"wp": workingpoint, "btag": "pass"}].view()
        ball = btag[{"wp": workingpoint, "btag": sum}].view()
        ball[ball<=0.]=1.
        ratio = bpass / np.maximum(ball, 1.)
        nom = hist.Hist(*btag.axes[2:], data=ratio)
        nom.name = "ratios"  
        nom.label = "out"
        self.eff = convert.from_histogram(nom).to_evaluator()

    def btag_weight(self, pt, eta, flavor, istag):

        abseta = abs(eta)
        flateta, counts = ak.fill_none(ak.flatten(abseta), 0.), ak.num(abseta)

        pt = ak.where((pt>999.99), ak.full_like(pt,999.99), pt)
        pt = ak.where((pt<20.0), ak.full_like(pt,20.0), pt)
        flatpt =  ak.fill_none(ak.flatten(pt), 20.)

        flatflavor = ak.fill_none(ak.flatten(flavor), 0)
        
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
        def P(eff):
            weight = ak.where(istag, eff, 1-eff)
            return ak.prod(weight, axis=1)

        '''
        particleNet_comb
         =>  systematic ,
         =>  working_point ,  L/M/T/XT/XXT
         =>  flavor ,  hadron flavor definition: 5=b, 4=c, 0=udsg
         =>  abseta ,
         =>  pt ,
        '''

        eff = ak.where(
            ~np.isnan(ak.fill_none(pt, np.nan)),
            ak.unflatten(self.eff.evaluate(flatflavor, flatpt, flateta), counts=counts),
            ak.zeros_like(pt)
        )
        print('wp', self._wp)
        print('flatflavor', ak.full_like(flatflavor, 0.))
        print('just flatflavor', flatflavor)
        print('flavor', flavor)
        print('flavor==0', (flavor==0)) 
        print('eta', flateta)
        print('pt' , flatpt)
        print('sf comb', self.sf['comb'])
        sf_nom = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_up_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_down_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_up_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_bc_down_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('central',self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('down_uncorrelated',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('down_uncorrelated',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)        
            )
        )
        sf_light_up_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('up_correlated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_down_correlated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('down_correlated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_up_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('up_uncorrelated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )
        sf_light_down_uncorrelated = ak.where(
            (flavor==0),
            ak.unflatten(self.sf['ligh'].evaluate('down_uncorrelated', self._wp, ak.full_like(flatflavor, 0.), flateta, flatpt), counts=counts),
            ak.where(
                (flavor==4),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 4.), flateta, flatpt), counts=counts),
                ak.unflatten(self.sf['comb'].evaluate('central',self._wp, ak.full_like(flatflavor, 5.), flateta, flatpt), counts=counts)
            )
        )

        eff_data_nom  = ak.where(
            (sf_nom*eff>1.), 
            ak.ones_like(eff), 
            sf_nom*eff
        )
        eff_data_bc_up_correlated   = ak.where(
            (sf_bc_up_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_up_correlated*eff
        )
        eff_data_bc_down_correlated = ak.where(
            (sf_bc_down_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_down_correlated*eff
        )
        eff_data_bc_up_uncorrelated = ak.where(
            (sf_bc_up_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_up_uncorrelated*eff
        )
        eff_data_bc_down_uncorrelated = ak.where(
            (sf_bc_down_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_bc_down_uncorrelated*eff
        )
        eff_data_light_up_correlated   = ak.where(
            (sf_light_up_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_up_correlated*eff
        )
        eff_data_light_down_correlated = ak.where(
            (sf_light_down_correlated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_down_correlated*eff
        )
        eff_data_light_up_uncorrelated = ak.where(
            (sf_light_up_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_up_uncorrelated*eff
        )
        eff_data_light_down_uncorrelated = ak.where(
            (sf_light_down_uncorrelated*eff>1.), 
            ak.ones_like(eff), 
            sf_light_down_uncorrelated*eff
        )

        nom = P(eff_data_nom)/P(eff)
        bc_up_correlated = P(eff_data_bc_up_correlated)/P(eff)
        bc_down_correlated = P(eff_data_bc_down_correlated)/P(eff)
        bc_up_uncorrelated = P(eff_data_bc_up_uncorrelated)/P(eff)
        bc_down_uncorrelated = P(eff_data_bc_down_uncorrelated)/P(eff)
        light_up_correlated = P(eff_data_light_up_correlated)/P(eff)
        light_down_correlated = P(eff_data_light_down_correlated)/P(eff)
        light_up_uncorrelated = P(eff_data_light_up_uncorrelated)/P(eff)
        light_down_uncorrelated = P(eff_data_light_down_uncorrelated)/P(eff)
        
        return np.nan_to_num(nom, nan=1.), \
        np.nan_to_num(bc_up_correlated, nan=1.), \
        np.nan_to_num(bc_down_correlated, nan=1.), \
        np.nan_to_num(bc_up_uncorrelated, nan=1.), \
        np.nan_to_num(bc_down_uncorrelated, nan=1.), \
        np.nan_to_num(light_up_correlated, nan=1.), \
        np.nan_to_num(light_down_correlated, nan=1.), \
        np.nan_to_num(light_up_uncorrelated, nan=1.), \
        np.nan_to_num(light_down_uncorrelated, nan=1.)

ids = load('data/ids.coffea')
corrections = load('data/corrections.coffea')
common = load('data/common.coffea')
isGoodAK4      = ids['isGoodAK4']
isGoodAK15      = ids['isGoodAK15']
isLooseElectron= ids['isLooseElectron']
isTightElectron= ids['isTightElectron']
isLooseMuon= ids['isLooseMuon']
isTightMuon= ids['isTightMuon']

PNetUParTWPs = common['btagWPs']['PNetUParT'][year]

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
mu_ntot = ak.num(mu, axis=1)
mu_nloose = ak.num(mu_loose, axis=1)
mu_ntight = ak.num(mu_tight, axis=1)
        

e = events.Electron
e['isloose'] = isLooseElectron(e,year)
e['istight'] = isTightElectron(e,year)
e['T'] = ak.zip(
    {
        "r": e.pt,
        "phi": e.phi,
    },
    with_name="PolarTwoVector",
    behavior=vector.behavior,
)

e_loose = e[e.isloose]
e_tight = e[e.istight]
e_ntot = ak.num(e, axis=1)
e_nloose = ak.num(e_loose, axis=1)
e_ntight = ak.num(e_tight, axis=1)

fj = events.AK15Puppi
fj['isgood'] = isGoodAK15(fj)
fj['T'] = ak.zip(
    {
        "r": fj.pt,
        "phi": fj.phi,
    },
    with_name="PolarTwoVector",
    behavior=vector.behavior,
)
fj_good = fj[fj.isgood]
fj_ntot = ak.num(fj, axis=1)
fj_ngood = ak.num(fj_good, axis=1)


j = events.Jet
j['T'] = ak.zip(
    {
        "r": j.pt,
        "phi": j.phi,
    },
    with_name="PolarTwoVector",
    behavior=vector.behavior,
)

j['isgood'] = isGoodAK4(j, year)
j['isbtagvL'] = (j.btagPNetB>PNetUParTWPs['loose'])

j_good = j[j.isgood]

j_ntot=ak.num(j, axis=1)
j_ngood=ak.num(j_good, axis=1)


print('j_good', j_good.pt, ak.type(j_good.pt))

btagSF, btagSFbc_correlatedUp, btagSFbc_correlatedDown, btagSFbc_uncorrelatedUp, btagSFbc_uncorrelatedDown, btagSFlight_correlatedUp, btagSFlight_correlatedDown, btagSFlight_uncorrelatedUp, btagSFlight_uncorrelatedDown  = BTagCorrector('PNetUParT',year,'loose').btag_weight( j_good.pt, j_good.eta, j_good.hadronFlavour, j_good.isbtagvL)
print('btagSF nom', btagSF)
print('btagSF bc correlated Up        :', btagSFbc_correlatedUp)
print('btagSF bc correlated Down      :', btagSFbc_correlatedDown)

print('btagSF bc uncorrelated Up      :', btagSFbc_uncorrelatedUp)
print('btagSF bc uncorrelated Down    :', btagSFbc_uncorrelatedDown)

print('btagSF light correlated Up     :', btagSFlight_correlatedUp)
print('btagSF light correlated Down   :', btagSFlight_correlatedDown)

print('btagSF light uncorrelated Up   :', btagSFlight_uncorrelatedUp)
print('btagSF light uncorrelated Down :', btagSFlight_uncorrelatedDown)
