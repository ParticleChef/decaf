import os, sys
import numpy as np
import matplotlib.pyplot as plt

fh = 'met_ah.txt'

nGen_fh = 44875968
nGen_sl = 54071479
nGen_dl = 56985050
lumi = 108950  # in pb-1
xsec_fh = 419.82 # in pb
xsec_sl = 405.75
xsec_dl = 98.04

data_fh = np.loadtxt(fh)
scale_fh = np.ones_like(data_fh) * lumi * xsec_fh / nGen_fh
data_sl = np.loadtxt('met_sl.txt')
scale_sl = np.ones_like(data_sl) * lumi * xsec_sl / nGen_sl
data_dl = np.loadtxt('met_test.txt')
scale_dl = np.ones_like(data_dl) * lumi * xsec_dl / nGen_dl
data = np.concatenate((data_fh, data_sl, data_dl))


bins = [250,260,270,280,290,300,350,400,500,800]
plt.style.use('seaborn-deep')
histo = plt.hist(data, bins=bins, weights= np.concatenate((scale_fh, scale_sl,  scale_dl)),
 histtype='stepfilled', alpha=1.0, color='yellow', edgecolor='black', label='TT')
print(histo[0])
plt.yscale('log')
plt.xlabel('Missing Transverse Energy (GeV)')
plt.ylabel('Events')
plt.ylim(0.1, 1e+8)
plt.grid(True, which='both', axis='y', ls='--', lw=0.5)
plt.legend( )
plt.savefig('f11.png')
