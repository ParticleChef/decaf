#! /usr/bin/env python
import correctionlib
import os
#import uproot_methods
import awkward as ak

import numpy as np
from coffea import hist, lookup_tools
from coffea.lookup_tools import extractor, dense_lookup

import uproot
from coffea import util
from coffea.util import save, load
import json




############
## BTag
## Btag recommend: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation#Recommendation_for_13_TeV_Data



from coffea.lookup_tools.correctionlib_wrapper import correctionlib_wrapper
from coffea.lookup_tools.dense_lookup import dense_lookup

class BTagCorrector:

    def __init__(self, tagger, year, workingpoint):
        self._year = year
        common = load('data/common.coffea')
        self._wp = common['btagWPs'][tagger][year][workingpoint]
        
        btvjson = correctionlib.CorrectionSet.from_file('data/BtagSF/'+year+'_UL/btagging.json.gz')
        self.sf = btvjson # ["deepJet_shape", "deepCSV_comb"]

        files = {
#            '2016preVFP': 'btageff2016.merged',
#            '2016postVFP': 'btageff2016.merged',
#            '2017': 'btageff2017.merged',
            '2018': 'btageff2018.merged',
        }
        filename = 'hists/'+files[year]
        btag = load(filename)
        bpass = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag', 'pass').values()[()]
        ball = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag').values()[()]
        ball[ball<=0.]=1.
        nom = bpass / np.maximum(ball, 1.)
        self.eff = lookup_tools.dense_lookup.dense_lookup(nom, [ax.edges() for ax in btag[tagger].axes()[3:]])

    
    def btag_weight(self, events, pt, eta, flavor, istag):
        tightJet = events
        abseta = abs(eta)
        
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1b_Event_reweighting_using_scale
        def P(eff):
            weight = np.ones_like(eff)
            weight_istag = ak.mask(eff, istag)
            weight_untag= ak.mask(1-eff, ~istag)
            weight = np.prod(weight_istag,axis=1)*np.prod(weight_untag,axis=1)

            #weight[~istag] = (1 - eff[~istag])
            return weight

        '''
        Correction deepJet_shape has 5 inputs
        Input systematic (string): 
        Input working_point (string): L/M/T
        Input flavor (int): hadron flavor definition: 5=b, 4=c, 0=udsg
        Input abseta (real):
        Input pt (real):
        '''
        #def bcSF(flav, pt, eta, syst):
        def bcSF(jetevents, flav, syst):
            Jet, nj = ak.flatten(jetevents), ak.num(jetevents)
            extsf = self.sf["deepJet_comb"].evaluate(syst, "M", np.ones_like(Jet.hadronFlavour) * flav , np.array(abs(Jet.eta)), np.array(Jet.pt))
            extsf = np.where(Jet.hadronFlavour == flav , extsf , 1.0)
            return ak.unflatten(extsf, nj)
    
        def lightSF(jetevents, flav, syst):
            Jet, nj = ak.flatten(jetevents), ak.num(jetevents)
            extsf = self.sf["deepJet_incl"].evaluate(syst, "M", np.ones_like(Jet.hadronFlavour) * flav , np.array(abs(Jet.eta)), np.array(Jet.pt))
            extsf = np.where(Jet.hadronFlavour == flav , extsf , 1.0)
            return ak.unflatten(extsf, nj)
    

        flav_b = flavor == 5
        flav_c = flavor == 4
        bc = flavor > 0
        light = ~bc
        
        eff = self.eff(flavor, pt, abseta)
        
        sf_b_nom = bcSF(tightJet, 5, 'central')
        sf_c_nom = bcSF(tightJet, 4, 'central')
        sf_l_nom = lightSF(tightJet, 0, 'central')
        #print('b', sf_b_nom)
        #print('c', sf_c_nom)
        #print('l', sf_l_nom)

        #eff_data_b_nom = np.minimum(1., sf_b_nom*eff)
        #eff_data_c_nom = np.minimum(1., sf_c_nom*eff)
        #eff_data_l_nom = np.minimum(1., sf_l_nom*eff)
        eff_data_b_nom = sf_b_nom*eff 
        eff_data_c_nom = sf_c_nom*eff
        eff_data_l_nom = sf_l_nom*eff



        nom = (P(eff_data_b_nom)*P(eff_data_c_nom)*P(eff_data_l_nom))/P(eff)

#        sf_nom = self.sf["deepJet_shape"].evaluate('central','M', flavor, abseta, pt)
        
#        bc_sf_up_correlated = pt.ones_like()
#        bc_sf_up_correlated[~bc] = sf_nom[~bc]
#        bc_sf_up_correlated[bc] = self.sf["deepJet_shape"].evaluate('up_correlated', 'M', flavor, eta, pt)[bc]
#        
#        bc_sf_down_correlated = pt.ones_like()
#        bc_sf_down_correlated[~bc] = sf_nom[~bc]
#        bc_sf_down_correlated[bc] = self.sf["deepJet_shape"].evaluate('down_correlated', 'M', flavor, eta, pt)[bc]
#
#        bc_sf_up_uncorrelated = pt.ones_like()
#        bc_sf_up_uncorrelated[~bc] = sf_nom[~bc]
#        bc_sf_up_uncorrelated[bc] = self.sf["deepJet_shape"].evaluate('up_uncorrelated', 'M', flavor, eta, pt)[bc]
#
#        bc_sf_down_uncorrelated = pt.ones_like()
#        bc_sf_down_uncorrelated[~bc] = sf_nom[~bc]
#        bc_sf_down_uncorrelated[bc] = self.sf["deepJet_shape"].evaluate('down_uncorrelated', 'M', flavor, eta, pt)[bc]
#
#        light_sf_up_correlated = pt.ones_like()
#        light_sf_up_correlated[~light] = sf_nom[~light]
#        light_sf_up_correlated[light] = self.sf["deepJet_shape"].evaluate('up_correlated', 'M', flavor, abseta, pt)[light]
#
#        light_sf_down_correlated = pt.ones_like()
#        light_sf_down_correlated[~light] = sf_nom[~light]
#        light_sf_down_correlated[light] = self.sf["deepJet_shape"].evaluate('down_correlated', 'M', flavor, abseta, pt)[light]
#
#        light_sf_up_uncorrelated = pt.ones_like()
#        light_sf_up_uncorrelated[~light] = sf_nom[~light]
#        light_sf_up_uncorrelated[light] = self.sf["deepJet_shape"].evaluate('up_uncorrelated', 'M', flavor, abseta, pt)[light]
#
#        light_sf_down_uncorrelated = pt.ones_like()
#        light_sf_down_uncorrelated[~light] = sf_nom[~light]
#        light_sf_down_uncorrelated[light] = self.sf["deepJet_shape"].evaluate('down_uncorrelated', 'M', flavor, abseta, pt)[light]



#        eff_data_nom  = np.minimum(1., sf_nom*eff)
#        bc_eff_data_up_correlated   = np.minimum(1., bc_sf_up_correlated*eff)
#        bc_eff_data_down_correlated = np.minimum(1., bc_sf_down_correlated*eff)
#        bc_eff_data_up_uncorrelated   = np.minimum(1., bc_sf_up_uncorrelated*eff)
#        bc_eff_data_down_uncorrelated = np.minimum(1., bc_sf_down_uncorrelated*eff)
#        light_eff_data_up_correlated   = np.minimum(1., light_sf_up_correlated*eff)
#        light_eff_data_down_correlated = np.minimum(1., light_sf_down_correlated*eff)
#        light_eff_data_up_uncorrelated   = np.minimum(1., light_sf_up_uncorrelated*eff)
#        light_eff_data_down_uncorrelated = np.minimum(1., light_sf_down_uncorrelated*eff)
       
#        nom = P(eff_data_nom)/P(eff)
#        bc_up_correlated = P(bc_eff_data_up_correlated)/P(eff)
#        bc_down_correlated = P(bc_eff_data_down_correlated)/P(eff)
#        bc_up_uncorrelated = P(bc_eff_data_up_uncorrelated)/P(eff)
#        bc_down_uncorrelated = P(bc_eff_data_down_uncorrelated)/P(eff)
#        light_up_correlated = P(light_eff_data_up_correlated)/P(eff)
#        light_down_correlated = P(light_eff_data_down_correlated)/P(eff)
#        light_up_uncorrelated = P(light_eff_data_up_uncorrelated)/P(eff)
#        light_down_uncorrelated = P(light_eff_data_down_uncorrelated)/P(eff)
        '''
        print('nom',sf_nom)
        print('bc_up_correlated',bc_sf_up_correlated)
        print('bc_down_correlated',bc_sf_down_correlated)
        print('bc_up_uncorrelated',bc_sf_up_uncorrelated)
        print('bc_down_uncorrelated',bc_sf_down_uncorrelated)
        print('light_up_correlated',light_sf_up_correlated)
        print('light_down_correlated',light_sf_down_correlated)
        print('light_up_uncorrelated',light_sf_up_uncorrelated)
        print('light_down_uncorrelated',light_sf_down_uncorrelated)
        '''
#        return np.nan_to_num(nom, nan=1.), np.nan_to_num(bc_up_correlated, nan=1.), np.nan_to_num(bc_down_correlated, nan=1.), np.nan_to_num(bc_up_uncorrelated, nan=1.), np.nan_to_num(bc_down_uncorrelated, nan=1.), np.nan_to_num(light_up_correlated, nan=1.), np.nan_to_num(light_down_correlated, nan=1.), np.nan_to_num(light_up_uncorrelated, nan=1.), np.nan_to_num(light_down_uncorrelated, nan=1.)
        return np.nan_to_num(nom, nan=1.)


get_btag_weight = {
    'deepflav': {
#        '2016preVFP': {
#            'loose'  : BTagCorrector('deepflav','2016preVFP','loose').btag_weight,
#            'medium' : BTagCorrector('deepflav','2016preVFP','medium').btag_weight,
#            'tight'  : BTagCorrector('deepflav','2016preVFP','tight').btag_weight
#        },
#        '2016postVFP': {
#            'loose'  : BTagCorrector('deepflav','2016postVFP','loose').btag_weight,
#            'medium' : BTagCorrector('deepflav','2016postVFP','medium').btag_weight,
#            'tight'  : BTagCorrector('deepflav','2016postVFP','tight').btag_weight
#        },
#        '2017': {
#            'loose'  : BTagCorrector('deepflav','2017','loose').btag_weight,
#            'medium' : BTagCorrector('deepflav','2017','medium').btag_weight,
#            'tight'  : BTagCorrector('deepflav','2017','tight').btag_weight
#        },
        '2018': {
            'loose'  : BTagCorrector('deepflav','2018','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2018','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2018','tight').btag_weight
        }
    },
    'deepcsv' : {
#        '2016preVFP': {
#            'loose'  : BTagCorrector('deepcsv','2016preVFP','loose').btag_weight,
#            'medium' : BTagCorrector('deepcsv','2016preVFP','medium').btag_weight,
#            'tight'  : BTagCorrector('deepcsv','2016preVFP','tight').btag_weight
#        },
#        '2016postVFP': {
#            'loose'  : BTagCorrector('deepcsv','2016postVFP','loose').btag_weight,
#            'medium' : BTagCorrector('deepcsv','2016postVFP','medium').btag_weight,
#            'tight'  : BTagCorrector('deepcsv','2016postVFP','tight').btag_weight
#        },
#        '2017': {
#            'loose'  : BTagCorrector('deepcsv','2017','loose').btag_weight,
#            'medium' : BTagCorrector('deepcsv','2017','medium').btag_weight,
#            'tight'  : BTagCorrector('deepcsv','2017','tight').btag_weight
#        },
        '2018': {
            'loose'  : BTagCorrector('deepcsv','2018','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2018','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2018','tight').btag_weight
        }
    }
}


class DoubleBTagCorrector:

    def __init__(self, year):
        self._year = year
        sf = {
            '2018': {
                'value': np.array([0.82, 0.82, 0.75, 0.81]),
                'unc': np.array([np.sqrt(0.07**2 + 0.11**2), np.sqrt(0.07**2 + 0.11**2), np.sqrt(0.06**2 + 0.06**2), np.sqrt(0.05**2 + 0.01**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
            '2017': {
                'value': np.array([0.84, 0.84, 0.98, 0.86]),
                'unc': np.array([np.sqrt(0.05**2 + 0.13**2), np.sqrt(0.05**2 + 0.13**2), np.sqrt(0.05**2 + 0.12**2), np.sqrt(0.05**2 + 0.05**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
            '2016': {
                'value': np.array([1.01, 1.01, 0.95, 0.99]),
                'unc': np.array([np.sqrt(0.06**2 + 0.02**2), np.sqrt(0.06**2 + 0.02**2), np.sqrt(0.05**2 + 0.09**2), np.sqrt(0.06**2 + 0.00**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
        }
        self.sf_nom=lookup_tools.dense_lookup.dense_lookup(sf[year]['value'], sf[year]['edges'])
        self.sf_up=lookup_tools.dense_lookup.dense_lookup(sf[year]['value']+sf[year]['unc'], sf[year]['edges'])
        self.sf_down=lookup_tools.dense_lookup.dense_lookup(sf[year]['value']-sf[year]['unc'], sf[year]['edges'])
        

    def doublebtag_weight(self, pt):
        return  self.sf_nom(pt), self.sf_up(pt), self.sf_down(pt)

get_doublebtag_weight = {
    '2016': DoubleBTagCorrector('2016').doublebtag_weight,
    '2017': DoubleBTagCorrector('2017').doublebtag_weight,
    '2018': DoubleBTagCorrector('2018').doublebtag_weight,
}



corrections = {}
corrections = {
    'get_btag_weight':          get_btag_weight,
}


save(corrections, 'data/correctionsBT.coffea')

