from coffea.util import save
import numpy as np

def match(a, b, val):
    combinations = a.cross(b, nested=True)
    return (combinations.i0.delta_r(combinations.i1)<val).any()

def sigmoid(x,a,b,c,d):
    """
    Sigmoid function for trigger turn-on fits.
    f(x) = c + (d-c) / (1 + np.exp(-a * (x-b)))
    """
    return c + (d-c) / (1 + np.exp(-a * (x-b)))


#########
## BTag wp
## UL18 WPs on ttbar RunIISummer19UL18MiniAOD dataset (pt>30 GeV)
## https://indico.cern.ch/event/967689/contributions/4083041/attachments/2130779/3590310/BTagPerf_201028_UL18WPs.pdf
## 2017 Ultra   Legacy  WPs on  RunIISummer19UL17MiniAOD    MC  production
## https://indico.cern.ch/event/880118/contributions/3714095/attachments/1973818/3284315/BTagPerf_200121_UL17WPs.pdf
## UL16 post-VFP WPs with RunIISummer20UL16MiniAOD* production (latest post-VPF)
## https://indico.cern.ch/event/1053254/contributions/4425737/attachments/2271561/3857890/BTagPerf_210623_UL16WPUpdate.pdf
## UL16 WPs separately for pre-VPF
## https://indico.cern.ch/event/1011636/contributions/4382495/attachments/2251885/3820785/BTagPerf_210525_UL16WPs.pdf

deepflavWPs = {
    '2016preVFP': {
        'loose' : 0.0508,
        'medium': 0.2598,
        'tight' : 0.6502
    },
    '2016postVFP': {
        'loose' : 0.0480,
        'medium': 0.2489,
        'tight' : 0.6377
    },
    '2017': {
        'loose' : 0.0532,
        'medium': 0.3040,
        'tight' : 0.7476
    },
    '2018': {
        'loose' : 0.0490,
        'medium': 0.2783,
        'tight' : 0.7100
    },
}
deepcsvWPs = {
    '2016preVFP': {
        'loose' : 0.2027,
        'medium': 0.6001,
        'tight' : 0.8819
    },
    '2016postVFP': {
        'loose' : 0.1918,
        'medium': 0.5847,
        'tight' : 0.8767
    },
    '2017': {
        'loose' : 0.1355,
        'medium': 0.4506,
        'tight' : 0.7738
    },
    '2018': {
        'loose' : 0.1208,
        'medium': 0.4168,
        'tight' : 0.7665
    },
}

btagWPs = {
    'deepflav': deepflavWPs,
    'deepcsv' : deepcsvWPs
}
    
common = {}
common['match'] = match
common['sigmoid'] = sigmoid
common['btagWPs'] = btagWPs
save(common, 'data/common.coffea')
