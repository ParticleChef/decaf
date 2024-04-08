#!/usr/bin/env python
'''@package docstring
Just a giant list of processes and properties
'''

processes =    {

    #data
    'MET':('MET','Data',1),
    'EGamma':('EGamma','Data',1),
    'SingleElectron':('SingleElectron','Data',1),
    'SinglePhoton':('SinglePhoton','Data',1),
    'SingleMuon':('SingleMuon','Data',1),
    
    # inclusive NLO V+jets 
#    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_nlo','MC',6025.2),
#    'DYJetsToNuNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZtoNuNu_nlo','MC',11433.),
#    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WJets_nlo','MC',61527.),
#    'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_lo_incl','MC',50400.),
#    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_lo_incl_CP5','MC',50400.),

    # LO Z->nunu
    #2018
    'Z1JetsToNuNu_M-50_LHEFilterPtZ-150To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8' :('Z1JetsToNuNu_M-50_LHEFilterPtZ-150To250', 'MC', 17.36),
    'Z1JetsToNuNu_M-50_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8' :('Z1JetsToNuNu_M-50_LHEFilterPtZ-250To400', 'MC', 1.978),
    'Z1JetsToNuNu_M-50_LHEFilterPtZ-400ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8' :('Z1JetsToNuNu_M-50_LHEFilterPtZ-400ToInf', 'MC', 0.2167),
    'Z1JetsToNuNu_M-50_LHEFilterPtZ-50To150_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'  :('Z1JetsToNuNu_M-50_LHEFilterPtZ-50To150', 'MC', 580.7),
    'Z2JetsToNuNu_M-50_LHEFilterPtZ-150To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8' :('Z2JetsToNuNu_M-50_LHEFilterPtZ-150To250', 'MC', 28.8),
    'Z2JetsToNuNu_M-50_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8' :('Z2JetsToNuNu_M-50_LHEFilterPtZ-250To400', 'MC', 4.989),
    'Z2JetsToNuNu_M-50_LHEFilterPtZ-400ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8' :('Z2JetsToNuNu_M-50_LHEFilterPtZ-400ToInf', 'MC', 0.8162),
    'Z2JetsToNuNu_M-50_LHEFilterPtZ-50To150_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'  :('Z2JetsToNuNu_M-50_LHEFilterPtZ-50To150', 'MC', 314.5),
    #2017
    #2016
    # Z->ll
    #2018
    'DYJetsToLL_LHEFilterPtZ-0To50_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'           :('DYJetsToLL_LHEFilterPtZ-0To50', 'MC', 1490.1),
    'DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'        :('DYJetsToLL_LHEFilterPtZ-100To250', 'MC', 94.39),
    'DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'        :('DYJetsToLL_LHEFilterPtZ-250To400', 'MC', 3.656),
    'DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'        :('DYJetsToLL_LHEFilterPtZ-400To650', 'MC', 0.4969),
    'DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'         :('DYJetsToLL_LHEFilterPtZ-50To100', 'MC', 395.1),
    'DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'        :('DYJetsToLL_LHEFilterPtZ-650ToInf', 'MC', 462.0),
    
    # W->lnu
    'WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8'                         :('WJetsToLNu_0J', 'MC', 53300.),
    'WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8'                         :('WJetsToLNu_1J', 'MC', 8947.0),
    'WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8'                         :('WJetsToLNu_2J', 'MC', 3335.0),
    'WJetsToLNu_Pt-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'   :('WJetsToLNu_Pt-100To250', 'MC', 757.7),
    'WJetsToLNu_Pt-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'   :('WJetsToLNu_Pt-250To400', 'MC', 27.53),
    'WJetsToLNu_Pt-400To600_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'   :('WJetsToLNu_Pt-400To600', 'MC', 3.511),
    'WJetsToLNu_Pt-600ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8'   :('WJetsToLNu_Pt-600ToInf', 'MC', 0.5426),

    #  gamma
    'G1Jet_LHEGpT-150To250_TuneCP5_13TeV-amcatnlo-pythia8'   :('G1Jet_LHEGpT-150To250', 'MC', 225.9),
    'G1Jet_LHEGpT-250To400_TuneCP5_13TeV-amcatnlo-pythia8'   :('G1Jet_LHEGpT-250To400', 'MC', 26.98),
    'G1Jet_LHEGpT-400To675_TuneCP5_13TeV-amcatnlo-pythia8'   :('G1Jet_LHEGpT-400To675', 'MC', 3.395),
    'G1Jet_LHEGpT-675ToInf_TuneCP5_13TeV-amcatnlo-pythia8'   :('G1Jet_LHEGpT-675ToInf', 'MC', 0.2477),
    
    # NLO W->lnu

    # NLO Z->ll

    # NLO Z->nunu

    # QCD
    'QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8'      :('QCD_Pt_1000to1400', 'MC', 7.482),
    'QCD_Pt_120to170_TuneCP5_13TeV_pythia8'        :('QCD_Pt_120to170', 'MC', 407500),
    'QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8'      :('QCD_Pt_1400to1800', 'MC', 0.1259),
    'QCD_Pt_15to30_TuneCP5_13TeV_pythia8'          :('QCD_Pt_15to30', 'MC', 1246000000.0),
    'QCD_Pt_170to300_TuneCP5_13TeV_pythia8'        :('QCD_Pt_170to300', 'MC', 103600),
    'QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8'      :('QCD_Pt_1800to2400', 'MC', 0.08748),
    'QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8'      :('QCD_Pt_2400to3200', 'MC', 0.005236),
    'QCD_Pt_300to470_TuneCP5_13TeV_pythia8'        :('QCD_Pt_300to470', 'MC', 6831),
    'QCD_Pt_30to50_TuneCP5_13TeV_pythia8'          :('QCD_Pt_30to50', 'MC', 106800000),
    'QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8'       :('QCD_Pt_3200toInf', 'MC', 0.0001351),
    'QCD_Pt_470to600_TuneCP5_13TeV_pythia8'        :('QCD_Pt_470to600', 'MC', 551.3),
    'QCD_Pt_50to80_TuneCP5_13TeV_pythia8'          :('QCD_Pt_50to80', 'MC', 15690000),
    'QCD_Pt_600to800_TuneCP5_13TeV_pythia8'        :('QCD_Pt_600to800', 'MC', 156.5),
    'QCD_Pt_800to1000_TuneCP5_13TeV_pythia8'       :('QCD_Pt_800to1000', 'MC', 26.15),
    'QCD_Pt_80to120_TuneCP5_13TeV_pythia8'         :('QCD_Pt_80to120', 'MC', 2341000),
    
    # Single tops
    'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8'                     :('ST_s-channel_4f_leptonDecays', 'MC', 6.96),
    'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8'    :('ST_t-channel_antitop_4f_InclusiveDecays', 'MC', 80.95),
    'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8'        :('ST_t-channel_top_4f_InclusiveDecays', 'MC', 136.02),
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8'                   :('ST_tW_antitop_5f_inclusiveDecays', 'MC', 35.85),
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8'                       :('ST_tW_top_5f_inclusiveDecays', 'MC', 35.85),


    # ttbar
    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8'           :('TTTo2L2Nu', 'MC', 88.3419),
    'TTToHadronic_TuneCP5_13TeV-powheg-pythia8'        :('TTToHadronic', 'MC', 377.9607),
    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'    :('TTToSemiLeptonic', 'MC', 365.4574),

        # regular dibosons
    'WW_TuneCP5_13TeV-pythia8'     :('WW_TuneCP5_13TeV-pythia8', 'MC', 119.),
    'WZ_TuneCP5_13TeV-pythia8'     :('WZ_TuneCP5_13TeV-pythia8', 'MC', 46.7),
    'ZZ_TuneCP5_13TeV-pythia8'     :('ZZ_TuneCP5_13TeV-pythia8', 'MC', 16.9),
    #'WW_TuneCP5_13TeV-pythia8':('Diboson_ww_CP5','MC',118.7),
    #'WZ_TuneCP5_13TeV-pythia8':('Diboson_wz_CP5','MC',47.13),
    #'ZZ_TuneCP5_13TeV-pythia8':('Diboson_zz_CP5','MC',16.523),


    # Private signal samples
    'TPhiTo2Chi_MPhi1000_MChi1000_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi1000_MChi1000', 'MC', 0.00005344),
    'TPhiTo2Chi_MPhi1000_MChi150_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi1000_MChi150', 'MC', 0.3727),
    'TPhiTo2Chi_MPhi1250_MChi150_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi1250_MChi150', 'MC', 0.1444),
    'TPhiTo2Chi_MPhi1495_MChi750_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi1495_MChi750', 'MC', 0.01206),
    'TPhiTo2Chi_MPhi1500_MChi1000_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi1500_MChi1000', 'MC', 0.0001071),
    'TPhiTo2Chi_MPhi1500_MChi150_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi1500_MChi150', 'MC', 0.06213),
    'TPhiTo2Chi_MPhi1700_MChi800_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi1700_MChi800', 'MC', 0.0189),
    'TPhiTo2Chi_MPhi1750_MChi150_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi1750_MChi150', 'MC', 0.02879.),
    'TPhiTo2Chi_MPhi1750_MChi700_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi1750_MChi700', 'MC', 0.02263),
    'TPhiTo2Chi_MPhi195_MChi100_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi195_MChi100', 'MC', 4.931),
    'TPhiTo2Chi_MPhi1995_MChi1000_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi1995_MChi1000', 'MC', 0.002655),
    'TPhiTo2Chi_MPhi2000_MChi1500_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi2000_MChi1500', 'MC', 0.000004319),
    'TPhiTo2Chi_MPhi2000_MChi150_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi2000_MChi150', 'MC', 0.01414),
    'TPhiTo2Chi_MPhi2000_MChi500_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi2000_MChi500', 'MC', 0.01297),
    'TPhiTo2Chi_MPhi200_MChi150_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi200_MChi150', 'MC', 0.2533),
    'TPhiTo2Chi_MPhi200_MChi50_TuneCP5_13TeV-amcatnlo-pythia8'   :('TPhiTo2Chi_MPhi200_MChi50', 'MC', 57.59),
    'TPhiTo2Chi_MPhi2495_MChi1250_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi2495_MChi1250', 'MC', 0.0006786),
    'TPhiTo2Chi_MPhi2500_MChi2000_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi2500_MChi2000', 'MC', 0.0000002682),
    'TPhiTo2Chi_MPhi2500_MChi750_TuneCP5_13TeV-amcatnlo-pythia8' :('TPhiTo2Chi_MPhi2500_MChi750', 'MC', 0.003237),
    'TPhiTo2Chi_MPhi295_MChi150_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi295_MChi150', 'MC', 2.373),
    'TPhiTo2Chi_MPhi2995_MChi1500_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi2995_MChi1500', 'MC', 0.0001901),
    'TPhiTo2Chi_MPhi3000_MChi1000_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi3000_MChi1000', 'MC', 0.0008935),
    'TPhiTo2Chi_MPhi3000_MChi2000_TuneCP5_13TeV-amcatnlo-pythia8':('TPhiTo2Chi_MPhi3000_MChi2000', 'MC', 0.0000004319),
    'TPhiTo2Chi_MPhi300_MChi100_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi300_MChi100', 'MC', 18.96),
    'TPhiTo2Chi_MPhi300_MChi300_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi300_MChi300', 'MC', 0.0189),
    'TPhiTo2Chi_MPhi495_MChi250_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi495_MChi250', 'MC', 0.704),
    'TPhiTo2Chi_MPhi500_MChi150_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi500_MChi150', 'MC', 4.315),
    'TPhiTo2Chi_MPhi500_MChi500_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi500_MChi500', 'MC', 0.002262),
    'TPhiTo2Chi_MPhi750_MChi150_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi750_MChi150', 'MC', 1.124),
    'TPhiTo2Chi_MPhi995_MChi500_TuneCP5_13TeV-amcatnlo-pythia8'  :('TPhiTo2Chi_MPhi995_MChi500', 'MC', 0.07102),



}
