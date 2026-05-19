#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from coffea.util import load
import mplhep
import matplotlib.pyplot as plt
from scipy.special import betaincinv  # type: ignore
import hist

myhist1 = load('hists/stop_2024.scaled')
myhist2 = load('hists/stop_2024_backup_251218.scaled')

bkg1 = myhist1['bkg']
bkg2 = myhist2['bkg']

print(
    'W (lnu)',
    bkg1['metpt']['W (lnu)'][{
        'region': 'cat1_preselection',
        'systematic': 'nominal'
    }].values()[-1]
)

print(
    'W (lnu)',
    bkg2['metpt']['W (lnu)'][{
        'region': 'cat1_preselection',
        'systematic': 'nominal'
    }].values()[-1]
)

def plot_compare(hist1, hist2, region, systematic, signal_name, bin_edges, xlabel, ylabel, title, filename):
    values1 = hist1['metpt'][signal_name][{
        'region': region,
        'systematic': systematic
    }].values()[-1]
    
    values2 = hist2['metpt'][signal_name][{
        'region': region,
        'systematic': systematic
    }].values()[-1]
    values1 = np.append(values1,values1[-1])
    values2 = np.append(values2,values2[-1])
    plt.style.use(mplhep.style.CMS)
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True, figsize=(9, 9))
    # Main plot
    ax.step(bin_edges, values2, where='post', label='WtoLNu Pt Bin', color='red', linewidth=2.5)
    ax.step(bin_edges, values1, where='post', label='WtoLNu Jet Bin', color='blue', linewidth=2.5)
    #h1 = plt.step(bin_edges, np.append(values2,values2[-1]), where='post', label='WtoLNu Pt Bin', color='blue', linewidth=2.5)
    #h2 = plt.step(bin_edges, np.append(values1,values1[-1]) , where='post', label='WtoLNu Jet Bin', color='red', linewidth=2.5)
    
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.grid()
    ax.set_xlim(250,800)
    ax.set_yscale('log')

    # Ratio plot
    ratio = np.divide(values2, values1, out=np.zeros_like(values1), where=values2!=0)
    rax.step(bin_edges, ratio, where='post', color='black', linewidth=2.5)
    rax.set_xlabel(xlabel)
    rax.set_ylabel('Pt/Jet')
    rax.grid()

    
    plt.savefig('wyield.png')
    plt.close()

bin_edges = bkg1['metpt']['W (lnu)'][{'region': 'cat1_preselection', 'systematic': 'nominal'}].axes[-1].edges
print(bin_edges)

plot_compare(
    bkg1,
    bkg2,
    region='cat1_preselection',
    systematic='nominal',
    signal_name='W (lnu)',
    bin_edges=bin_edges,
    xlabel='MET pT (GeV)',
    ylabel='Events',
    title='Comparison of W (lnu) Yields',
    filename='wyield.png'
)
