import hist
from coffea.util import load, save
import mplhep as hep
import matplotlib.pyplot as plt
import copy
import numpy as np
import os
from tabulate import tabulate
import uncertainties as unc  
import uncertainties.unumpy as unumpy

year = "2022post"
histversion = year #'2022pre'#22Test'#'GenIsoPhoTest' # 'sdfjTest'#
whathist = {
    '2022pre': '2022_0616_test1',
    '2022post':'2022EE_0616',
    '2023pre': '2023_0616',
    '2023post': '2023BPix_0616'

}
lumi = {
    '2018': 59.83,
    '2022pre' : 7.99,
    '2022post': 26.68,
    '2023pre': 17.96,
    '2023post': 9.68
}
if 'Test' in histversion:
    hists = load('./hadmonotop'+whathist[histversion]+'.scaled')
else:
    hists = load('./hadmonotop'+whathist[histversion]+'.scaled')
    #hists = load('./hadmonotop'+year+'_'+whathist[histversion]+'.scaled')

hists_bkg = hists['bkg']
hists_data= hists['data']


############
## Yields
############
TPhi = '' #'TPhiTo2Chi_MPhi200_MChi150_TuneCP5_13TeV-amcatnlo-pythia8'
variables = list(hists_bkg.keys())
process = list(hists_bkg['sumw'].keys())
print(process)
region = 'wmcr'
regions = ['sr','tecr', 'tmcr', 'wecr', 'wmcr', 'zecr', 'zmcr', 'gcr']
print(f'** Yields table for {year} {lumi[year]} fb-1 **')
v = 'ut'
for reg in regions:

    if 's' in reg:
        data_name = 'MET'
    elif 'e' in reg:
        data_name = 'EGamma'
        #continue
    elif 'm' in reg:
        data_name = 'MET'
    else:
        data_name = 'EGamma'
        #continue
    bkg_all = 0
    print('======== ', reg, ' ========')
    for proc in process:
        #if 'Z' in proc : continue
        if reg not in list(hists_bkg[v][proc].axes[0]) :
            print("%-25s" % proc, "%10.2f" % 0.0)
        else:
            print("%-25s" % proc, "%10.2f" % np.sum(hists_bkg[v][proc][{'region': reg}].values()))
            bkg_all = bkg_all + np.sum(hists_bkg[v][proc][{'region': reg}].values())
    print('------------------------------------')
    print("%-25s" % str('BKG '), "%10.2f" % bkg_all)
    print("%-25s" % str('DATA'), "%10.2f" % np.sum(hists_data[v][data_name][{'region': reg}].values()))
    print()
