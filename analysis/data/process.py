#!/usr/bin/env python

processes = {

    ## Data
    'JetMET':('JetMET','Data',1),
    'EGamma':('EGamma','Data',1),

    ## DYto2L-2Jets
    'DYto2L-2Jets_MLL-50_PTLL-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-100to200_1J', 'MC', '45.42'),
    'DYto2L-2Jets_MLL-50_PTLL-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-100to200_2J', 'MC', '51.68'),
    'DYto2L-2Jets_MLL-50_PTLL-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-200to400_1J', 'MC', '3.382'),
    'DYto2L-2Jets_MLL-50_PTLL-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-200to400_2J', 'MC', '7.159'),
    'DYto2L-2Jets_MLL-50_PTLL-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-400to600_1J', 'MC', '0.1162'),
    'DYto2L-2Jets_MLL-50_PTLL-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-400to600_2J', 'MC', '0.4157'),
    'DYto2L-2Jets_MLL-50_PTLL-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-40to100_1J', 'MC', '475.3'),
    'DYto2L-2Jets_MLL-50_PTLL-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-40to100_2J', 'MC', '179.3'),
    'DYto2L-2Jets_MLL-50_PTLL-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-600_1J', 'MC', '0.01392'),
    'DYto2L-2Jets_MLL-50_PTLL-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('DYto2L-2Jets_MLL-50_PTLL-600_2J', 'MC', '0.07019'),

    ## DYto2L-4Jets
    'DYto2L-4Jets_MLL-50_PTLL-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8': ('DYto2L-4Jets_MLL-50_PTLL-100to200', 'MC', '58.46'),
    'DYto2L-4Jets_MLL-50_PTLL-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8': ('DYto2L-4Jets_MLL-50_PTLL-200to400', 'MC', '6.678'),
    'DYto2L-4Jets_MLL-50_PTLL-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8': ('DYto2L-4Jets_MLL-50_PTLL-400to600', 'MC', '0.3833'),
    'DYto2L-4Jets_MLL-50_PTLL-40to100_TuneCP5_13p6TeV_madgraphMLM-pythia8': ('DYto2L-4Jets_MLL-50_PTLL-40to100', 'MC', '403.7'),
    'DYto2L-4Jets_MLL-50_PTLL-600_TuneCP5_13p6TeV_madgraphMLM-pythia8': ('DYto2L-4Jets_MLL-50_PTLL-600', 'MC', '0.06843'),

    ## GJ
    'GJ_PTG-100to200_TuneCP5_13p6TeV_amcatnlo-pythia8': ('GJ_PTG-100to200', 'MC', '1396.0'),
    'GJ_PTG-200to400_TuneCP5_13p6TeV_amcatnlo-pythia8': ('GJ_PTG-200to400', 'MC', '88.52'),
    'GJ_PTG-400to600_TuneCP5_13p6TeV_amcatnlo-pythia8': ('GJ_PTG-400to600', 'MC', '3.783'),
    'GJ_PTG-600_TuneCP5_13p6TeV_amcatnlo-pythia8': ('GJ_PTG-600', 'MC', '0.5755'),

    ## QCD
    'QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8': ('QCD_PT-1000to1400', 'MC', '8.92'),
    'QCD_PT-120to170_TuneCP5_13p6TeV_pythia8': ('QCD_PT-120to170', 'MC', '445800.0'),
    'QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8': ('QCD_PT-1400to1800', 'MC', '0.8103'),
    'QCD_PT-15to30_TuneCP5_13p6TeV_pythia8': ('QCD_PT-15to30', 'MC', '1301000000.0'),
    'QCD_PT-170to300_TuneCP5_13p6TeV_pythia8': ('QCD_PT-170to300', 'MC', '113700.0'),
    'QCD_PT-300to470_TuneCP5_13p6TeV_pythia8': ('QCD_PT-300to470', 'MC', '7559.0'),
    'QCD_PT-30to50_TuneCP5_13p6TeV_pythia8': ('QCD_PT-30to50', 'MC', '113300000.0'),
    'QCD_PT-470to600_TuneCP5_13p6TeV_pythia8': ('QCD_PT-470to600', 'MC', '626.4'),
    'QCD_PT-50to80_TuneCP5_13p6TeV_pythia8': ('QCD_PT-50to80', 'MC', '16760000.0'),
    'QCD_PT-600to800_TuneCP5_13p6TeV_pythia8': ('QCD_PT-600to800', 'MC', '178.6'),
    'QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8': ('QCD_PT-800to1000', 'MC', '30.57'),
    'QCD_PT-80to120_TuneCP5_13p6TeV_pythia8': ('QCD_PT-80to120', 'MC', '2534000.0'),

    ## TbarBtoLminusNuB-s-channel-4FS
    'TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8': ('TbarBtoLminusNuB-s-channel-4FS', 'MC', '1.43'),

    ## TbarWplusto2L2Nu
    'TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8': ('TbarWplusto2L2Nu', 'MC', '3.9735'),

    ## TbarWplusto4Q
    'TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8': ('TbarWplusto4Q', 'MC', '15.942'),

    ## TbarWplustoLNu2Q
    'TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8': ('TbarWplustoLNu2Q', 'MC', '15.918'),

    ## TBbartoLplusNuBbar-s-channel-4FS
    'TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8': ('TBbartoLplusNuBbar-s-channel-4FS', 'MC', '2.278'),

    ## TQbarto2Q-t-channel
    'TQbarto2Q-t-channel_TuneCP5_13p6TeV_powheg-pythia8': ('TQbarto2Q-t-channel', 'MC', '43.16'),

    ## TQbartoLNu-t-channel
    'TQbartoLNu-t-channel_TuneCP5_13p6TeV_powheg-pythia8': ('TQbartoLNu-t-channel', 'MC', '86.45'),

    ## TTto2L2Nu
    'TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8': ('TTto2L2Nu', 'MC', '101.802'),

    ## TTto4Q
    'TTto4Q_TuneCP5_13p6TeV_powheg-pythia8': ('TTto4Q', 'MC', '408.439'),

    ## TTtoLNu2Q
    'TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8': ('TTtoLNu2Q', 'MC', '407.824'),

    ## TWminusto2L2Nu
    'TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8': ('TWminusto2L2Nu', 'MC', '3.966'),

    ## TWminusto4Q
    'TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8': ('TWminusto4Q', 'MC', '15.915'),

    ## TWminustoLNu2Q
    'TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8': ('TWminustoLNu2Q', 'MC', '15.891'),

    ## WtoLNu-2Jets
    'WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_0J', 'MC', '55760.0'),
    'WtoLNu-2Jets_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_1J', 'MC', '9529.0'),
    'WtoLNu-2Jets_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_2J', 'MC', '3532.0'),
    'WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-100to200_1J', 'MC', '368.2'),
    'WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-100to200_2J', 'MC', '421.9'),
    'WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-200to400_1J', 'MC', '25.6'),
    'WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-200to400_2J', 'MC', '54.77'),
    'WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-400to600_1J', 'MC', '0.8785'),
    'WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-400to600_2J', 'MC', '3.119'),
    'WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-40to100_1J', 'MC', '4427.0'),
    'WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-40to100_2J', 'MC', '1598.0'),
    'WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-600_1J', 'MC', '0.1053'),
    'WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('WtoLNu-2Jets_PTLNu-600_2J', 'MC', '0.5261'),

    ## WW
    'WW_TuneCP5_13p6TeV_pythia8': ('WW', 'MC', '80.23'),

    ## WZ
    'WZ_TuneCP5_13p6TeV_pythia8': ('WZ', 'MC', '29.1'),

    ## Zto2Nu-2Jets
    'Zto2Nu-2Jets_PTNuNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-100to200_1J', 'MC', '86.38'),
    'Zto2Nu-2Jets_PTNuNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-100to200_2J', 'MC', '100.4'),
    'Zto2Nu-2Jets_PTNuNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-200to400_1J', 'MC', '6.354'),
    'Zto2Nu-2Jets_PTNuNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-200to400_2J', 'MC', '13.86'),
    'Zto2Nu-2Jets_PTNuNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-400to600_1J', 'MC', '0.2188'),
    'Zto2Nu-2Jets_PTNuNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-400to600_2J', 'MC', '0.7816'),
    'Zto2Nu-2Jets_PTNuNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-40to100_1J', 'MC', '929.8'),
    'Zto2Nu-2Jets_PTNuNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-40to100_2J', 'MC', '335.5'),
    'Zto2Nu-2Jets_PTNuNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-600_1J', 'MC', '0.02583'),
    'Zto2Nu-2Jets_PTNuNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8': ('Zto2Nu-2Jets_PTNuNu-600_2J', 'MC', '0.1311'),

    ## ZZ
    'ZZ_TuneCP5_13p6TeV_pythia8': ('ZZ', 'MC', '12.75'),
}
