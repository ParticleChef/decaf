import uproot
import numpy as np
import awkward as ak
from coffea.util import load, save
import mplhep
import matplotlib.pyplot as plt
import hist
from scipy.special import betaincinv # type: ignore
error_opts = {
    'step': 'post',
    'label': 'Stat. unc',
    'hatch': '///',
    'facecolor': 'black',
    'alpha': 0.3,
    'edgecolor': (0, 0, 0, 0.3),
    'linewidth': 0
}

myhist = load('hists/stop_2024.scaled')
bkg = myhist['bkg']
data = myhist['data']

print('TT',bkg['metpt']['TT'][{'region':'cat1_preselection','systematic':'nominal'}].values())
#print('Z (inv)',bkg['metpt']['Z (inv)'][{'region':'cat1_preselection','systematic':'nominal'}].values())
#print('W (lnu)',bkg['metpt']['W (lnu)'][{'region':'cat1_preselection','systematic':'nominal'}].values())
print('Single Top',bkg['metpt']['Single Top'][{'region':'cat1_preselection','systematic':'nominal'}].values())
print('QCD Multijet',bkg['metpt']['QCD Multijet'][{'region':'cat1_preselection','systematic':'nominal'}].values())
print('data',data['metpt']['JetMET'][{'region':'cat1_preselection','systematic':'nominal'}].values())

regions = 'cat1_preselection'
systematics = 'nominal'
#stacks = ['QCD Multijet','Z (inv)','W (lnu)','Single Top', 'TT']
#colors= ['lightblue','salmon','green','orange','yellow']
#stacks = ['Z (inv)', 'W (lnu)', 'TT','Single Top','QCD Multijet']
#colors = ['salmon','green','yellow','orange','lightblue']
stacks = ['QCD Multijet','TT', 'Single Top']
colors = ['lightblue', 'yellow', 'orange']


for key in bkg.keys():
    if 'sumw' in key:
        continue
    elif 'template' in key: continue
    #if 'n' in key: continue
    else:
        print("Plotting ",key)
        plt.style.use(mplhep.style.CMS)
        fig, (ax, rax) = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(9, 9),
            gridspec_kw={"height_ratios": (3, 1)},
            sharex=True,
        )
        fig.subplots_adjust(hspace=0.07)
        mplhep.cms.label(ax=ax, llabel='Work in progress', rlabel='108.95 fb$^{-1}$ (13.6 TeV)', fontsize=24)
        ax.set_prop_cycle(color=colors)
        #print("drawing ",key)
        for mc in stacks:
            # stacked histograms
            h = bkg[key][mc][{'region': regions, 'systematic': systematics}]
            try:
                # just for metpt
                h = h[{'signal_trigger': 0}] + h[{'signal_trigger': 1}]
            except:
                pass
            # stacking
            c = h.copy()
            try:
                origin_stack += c
                #origin_stack += h
                draw_stack += h.values()
                
            except:
                origin_stack = c
                #origin_stack = h
                draw_stack = h.values()
            #print(origin_stack.view().value)
        bins = h.axes[0].edges
        #print(bins)

        total_bkg = np.array(origin_stack.view().value)
        unc_low = total_bkg - np.sqrt(origin_stack.view().variance)
        unc_up = total_bkg + np.sqrt(origin_stack.view().variance)
        ## append dummy value for last step
        total_bkg_forratio = np.append(total_bkg, total_bkg[-1])
        unc_low = np.append(unc_low, unc_low[-1])
        unc_up = np.append(unc_up, unc_up[-1])

        for mc in stacks:
            #print(mc, draw_stack)
            # stacked histograms
            if 'Z (inv)' in mc:
                labels = r'Z$\rightarrow$$\nu$$\bar{\nu}$'
            elif 'W (lnu)' in mc:
                labels = r'W$\rightarrow$$\ell$$\nu$'
            elif 'TT' in mc:
                labels = 't$\\bar{t}$'
            #elif 'ST' in mc:
            #    labels = 'Single Top'
            else:
                labels = mc
            #print(len(bins[:-1]), len(draw_stack))
            ax.hist(bins[:-1], bins, weights=draw_stack, histtype='stepfilled', label=labels, alpha=1.0, edgecolor='black')
            temp = bkg[key][mc][{'region': regions, 'systematic': systematics}]
            try:
                # just for metpt
                temp = temp[{'signal_trigger': 0}] + temp[{'signal_trigger': 1}]
                #print(mc, temp.values())
            except:
                pass
            draw_stack = draw_stack - temp.values()

        ## uncertainty band
        
        ax.fill_between(bins, unc_low, unc_up, **error_opts)
        
        # data points with error bars
        h_data = data[key]['JetMET'][{'region': regions, 'systematic': systematics}]
        try:
            h_data = h_data[{'signal_trigger': 0}] + h_data[{'signal_trigger': 1}]
        except:
            pass
        bins = h_data.axes[0].edges
        centers = h_data.axes[0].centers
        ax.errorbar(
            centers, h_data.values(),
            xerr = np.diff(bins) / 2,
            yerr = np.sqrt(h_data.values()),
            fmt='o', label='Data',
            markersize=8, capsize=5, capthick=1,
            color='black',
        )
        ax.set_ylabel('Events')
        ax.legend(ncol=2)
        ax.set_yscale('log')
        ax.set_yticks([0.1, 1, 10, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8])
        ax.set_ylim(0.1, 1e+8)
        if 'metpt' in key:
            ax.set_xlim(250, 800)
        if 'eta' in key:
            ax.set_xlim(-3.0, 3.0)
        if 'phi' in key:
            ax.set_xlim(-3.2, 3.2)
        if 'nElectron' in key or 'nMuon' in key:
            ax.set_xlim(0,4)
        if 'nJet' in key:
            ax.set_xlim(0,10)
        ax.grid(True, which='both', axis='y', ls='--', lw=0.5)
        
        # ratio plot
        rax.errorbar(
            centers, h_data.values()/total_bkg,
            xerr = np.diff(bins) / 2,
            yerr = np.sqrt(h_data.values())/total_bkg,
            fmt='o',
            color='black',
            markersize=8, capsize=5, capthick=1,
        )
        #print(unc_low, unc_up)
        rax.fill_between(bins, (unc_low/total_bkg_forratio), (unc_up/total_bkg_forratio), **error_opts)
        rax.set_ylabel('Data/MC')
        rax.set_ylim(0, 2)
        rax.set_yticks([0, 0.5, 1, 1.5, 2])
        rax.grid(True, which='both', axis='y', ls='--', lw=0.5)
        plt.xlabel(key)
        plt.savefig('plots/'+key+'_'+regions+'.png')
        plt.close()

