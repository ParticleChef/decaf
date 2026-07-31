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


year = "2022pre"
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



regions = ['sr', 'tecr', 'tmcr', 'wecr', 'wmcr', 'zecr', 'zmcr', 'gcr']
dataset = 'Z ($\\nu\\nu$) + Jets'
data_cutflow = hists_bkg['cutflow']
print(f'** {dataset} cutflow for {year} {lumi[year]} fb-1 **')
for reg in regions:

    cutlists=data_cutflow[dataset].axes[1]
    yd = hists['bkg']['cutflow'][dataset][{'region':reg}].values()[()]
    real_yd = np.sum(yd, axis=1)
    show_me = {}
    for idx, cut in enumerate(cutlists):
        if real_yd[idx] == 0.0:
            #print('value zero %s' % cut)
            continue
        else:
            show_me[cut] = real_yd[idx]


    real_cutflow = {k: v for k, v in sorted(show_me.items(), key=lambda item: item[1], reverse=True)}

    print()
    #print("== Hist version is ", whathist[histversion], " ==")
    print("          == [", reg, "] ==")
    for key in real_cutflow:
        if 'exclude' in key: continue
        print("%-27s" % key, "%13.2f" % real_cutflow[key])
    print()
