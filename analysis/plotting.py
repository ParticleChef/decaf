import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import coffea
from coffea import util
import hist
import os, sys
import argparse
from tabulate import tabulate
import uncertainties as unc  
import uncertainties.unumpy as unp
from cycler import cycler
## read merged file
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="input file")
parser.add_argument("-y", type=str, help="year", default="2018", choices=["2016", "2017", "2018"])
parser.add_argument("-b", type=str, help="datablind", default="", choices=["", "blind"])
args = parser.parse_args()

class Mydrawer:
    def __init__(self, f, y, b):
        try:
            self.f = util.load(f)
            self.y = y
            self.b = b
        except:
            print("Please provide a valid file. Exiting...")
            sys.exit(1)
        self.fname = str(f).split("/")[-1].split(".")[0]
        self.ftype = str(f).split("/")[-1].split(".")[1]
        self.workingdir = os.getcwd()
        self.outdir = self.workingdir + "/plots/" + self.fname
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
    
    def scaleDraw(self):
        lumis = {'2016': 35.9, '2017': 41.5, '2018': 59.8}
        ftype = self.ftype
        if ftype != "scaled":
            print("This is not a scaled file. Exiting...")
            sys.exit(1)
        hists = self.f
        year = self.y
        blind = self.b
        bkg = hists['bkg']
        sig = hists['sig']
        data = hists['data']
        regions = ['sr', 'wmcr', 'wecr', 'tmcr', 'gcr', 'tecr', 'zecr', 'zmcr']
        stacks = {
            'sr':[r'Z ($\nu\nu$) + Jets','TT',r'W ($\ell\nu$) + Jets', 'QCD Multijet', 'VV','Single Top', 'Z ($\ell\ell$) + Jets', '$\gamma$ + Jets'],
            'tw':[r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV','Single Top', 'Z ($\ell\ell$) + Jets',r'Z ($\nu\nu$) + Jets', '$\gamma$ + Jets'],
            'zr':['Z ($\ell\ell$) + Jets',r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV','Single Top', r'Z ($\nu\nu$) + Jets', '$\gamma$ + Jets'],
            'gr':['$\gamma$ + Jets',r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV','Single Top', 'Z ($\ell\ell$) + Jets',r'Z ($\nu\nu$) + Jets'],
        }
        #stacks = [r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV','Single Top', 'Z ($\ell\ell$) + Jets',r'Z ($\nu\nu$) + Jets']
        colors = {
            'sr':['gold','blue','maroon','green','lightgreen','orange','brown','purple'],
            'tw':['maroon','blue','green','lightgreen','orange','brown','gold','purple'],
            'zr':['brown','maroon','blue','green','lightgreen','orange','gold','purple'],
            'gr':['purple','blue','maroon','green','lightgreen','orange','brown','gold'],
        }
        leg_reg = {'sr': 'SR', 'wmcr': r'W($\mu$) CR', 'wecr': r'W(e) CR', 'tmcr': r'Top($\mu$) CR', 'gcr': 'G CR', 'tecr': r'Top(e) CR', 'zecr': r'Z(ee) CR', 'zmcr': r'Z($\mu\mu$) CR'}
        for region in regions:
            if 'sr' in region: stack_map = 'sr'
            elif 'z' in region: stack_map = 'zr'
            elif 'g' in region: stack_map = 'gr'
            else: stack_map = 'tw'
            print("-----------------------------------------------------")
            for key in bkg.keys():
                print(region, key)
                if len(bkg[key]) == 0:
                    print(region, key, 'is empty')
                    continue
                if 'sumw' in key:
                    continue
                plt.style.use(hep.style.CMS)
                
                fig, (ax, rax) = plt.subplots(
                    nrows=2,
                    ncols=1,
                    figsize=(10,10),
                    gridspec_kw={"height_ratios": (3, 1)},
                    sharex=True
                )
                fig.subplots_adjust(hspace=.07)
                hep.cms.label(ax=ax, llabel='Private Work', rlabel=str(lumis[year])+' fb$^{-1}$ (13 TeV)')
                ax.set_prop_cycle(cycler(color=colors[stack_map]))
                
                ### Plotting loop
                for sample in stacks[stack_map]:
                    try:
                        target = bkg[key][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}]

                        bins = target.axes.edges[0]
                        try:
                            mcstack += bkg[key][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}]
                            drawstack += bkg[key][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].values()
                        except:
                            mcstack = bkg[key][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}]
                            drawstack = bkg[key][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].values()
                    except:
                        continue                
                
                if region == 'sr':
                    sig_key = 'TPhiTo2Chi_MPhi1000_MChi150_TuneCP5_13TeV-amcatnlo-pythia8'
                    try:
                        sig[key][sig_key][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].plot1d(ax=ax, histtype='step', color='c', stack=False, label='Monotop_MPhi1500_MChi150')
                        sigstack = sig[key][sig_key][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}]
                    except:
                        print("Not able to draw... skipping", key, 'Data')
                        continue

                for sample in stacks[stack_map]:
                    try:
                        ax.hist(bins[:-1], bins, weights=drawstack, histtype='stepfilled', label=sample, linewidth=2, edgecolor=(0,0,0,0.3), alpha=1.0)
                        drawstack -= bkg[key][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].values()
                    except:
                        continue

                ### Uncertainty histogram as bin by bin
                error_opts = {
                    'step': 'post',
                    'label': 'Stat. Unc.',
                    'hatch': '///',
                    'facecolor': 'None',
                    'edgecolor': (0, 0, 0, 0.3),
                    'linewidth': 0
                }
                try:
                    unc_low =  mcstack.view().value - np.sqrt(mcstack.view().variance)
                    unc_high = mcstack.view().value + np.sqrt(mcstack.view().variance)
                    ax.fill_between(bins[:-1], unc_low, unc_high, **error_opts)
                except:
                    pass

                if blind != 'blind':
                    if region == 'sr' or 'mcr' in region:
                        try:
                            data[key]['MET'][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].plot1d(ax=ax, histtype='errorbar', color='k', stack=False, label='Data')
                            datastack = data[key]['MET'][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}]
                        except:
                            print("Not able to draw... skipping", key, 'Data')
                            continue

                    elif region == 'gcr' or 'ecr' in region:
                        try:
                            data[key]['EGamma'][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].plot1d(ax=ax, histtype='errorbar', color='k', stack=False, label='Data')
                            datastack = data[key]['EGamma'][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}]
                        except:
                            print("Not able to draw... skipping", key, 'Data')
                            continue

                ax.set_yscale('log')
                ax.set_ylim(0.01, 1000000)
                if key == 'met':
                    ax.set_xlim(0, 1200)
                    rax.set_xlabel(r'E$_{T}$ (GeV)')
                leg = ax.legend(ncol=3, loc='upper left', fontsize=16, title= leg_reg[region])

                def legend_title_left(leg):
                    c = leg.get_children()[0]
                    title = c.get_children()[0]
                    hpack = c.get_children()[1]
                    c._children = [hpack]
                    hpack._children = [title] + hpack.get_children()

                legend_title_left(leg)
                leg.get_title().set_fontsize('16')

                ax.set_xlabel(None)
                ax.set_ylabel('Events / bins')
                ax.grid(True)
                
                # Draw the ratio plot
                if blind != 'blind':
                    ds = datastack.view().value
                    ms = mcstack.view().value
                    ratio = ds / ms
                    rax.errorbar(
                        x=target.axes.edges[0][:-1] + np.diff(target.axes.edges[0]) / 2,
                        y=ratio,
                        yerr=np.sqrt(datastack.view().variance) / mcstack.view().value,
                        fmt="o",
                        color="k",
                    )
                    rax.fill_between(
                        bins,
                        # append the last edge to the bins
                        (1.0 - np.sqrt(mcstack.view().variance) / mcstack.view().value).tolist() + [1.0],
                        (1.0 + np.sqrt(mcstack.view().variance) / mcstack.view().value).tolist() + [1.0],
                        step="post",
                        color="k",
                        alpha=0.5,
                    )
                if not key == 'template':
                    lab = bkg[key][sample].axes[key].label
                    rax.set_xlabel(lab)
                rax.set_ylabel('Data/MC')
                rax.set_ylim(0.5, 1.5)
                rax.axhline(1, color='k', linestyle='--',alpha=0.5)
                rax.grid(True)
                fig.savefig(self.outdir+'/'+region+'_'+key+'.png')
                plt.close()

            # print yield using fj1pt
            print("-----------------------------------------------------")
            print(region, 'yield')
            sum_stack = 0
            for sample in stacks:
                try:
                    # print upper stat unc.
                    print("{0:<25} {1:>10}".format(sample, np.sum(bkg['fj1pt'][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].values())))
                    sum_stack += np.sum(bkg['fj1pt'][sample][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].values())
                except:
                    continue
                
            print("-----------------------------------------------------")
            print("{0:<25} {1:>10}".format('Total BKG', np.sum(mcstack.view().value)))
            if blind != 'blind':
                if region == 'sr' or 'mcr' in region:
                    print("{0:<25} {1:>10}".format('Data', np.sum(data['fj1pt']['MET'][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].values())))
                elif region == 'gcr' or 'ecr' in region:
                    print("{0:<25} {1:>10}".format('Data', np.sum(data['fj1pt']['EGamma'][{'TvsQCD': slice(0.0j,0.26j,sum), 'region': region}].values())))

    def btagSpliter(self, wp):
        hists = self.f
        year = self.y
        keys = hists.keys()
        deepflav = hists["deepflav"]
        deepcsv = hists["deepcsv"]
        if year == "2016":
            btag = deepcsv
        else:
            btag = deepflav
        for mcset in hists['deepflav']:
            try:
                btag += hists['deepflav'][mcset]
            except:
                btag = hists['deepflav'][mcset]
            bpass = btag[{"wp": str(wp), "btag": "pass"}].view()
            ball = btag[{"wp": str(wp), "btag": sum}].view()
            upass = bpass
            uall = ball
            # Light Flavor Jets
            lfpass = unp.uarray((bpass[0], np.sqrt(upass[0])))
            lfall = unp.uarray((ball[0], np.sqrt(uall[0])))
            lfall[lfall <= 0] = 1.00
            lf = lfpass / lfall
            # Charm Flavor Jets
            cfpass = unp.uarray((bpass[1], np.sqrt(upass[1])))
            cfall = unp.uarray((ball[1], np.sqrt(uall[1])))
            cfall[cfall <= 0] = 1.00
            cf = cfpass / cfall
            # Bottom Flavor Jets
            bfpass = unp.uarray((bpass[2], np.sqrt(upass[2])))
            bfall = unp.uarray((ball[2], np.sqrt(uall[2])))
            bfall[bfall <= 0] = 1.00
            bf = bfpass / bfall
            return lf, cf, bf

    def btagDraw(self, wp):
        year = self.y
        yedges = [30, 50, 70, 100, 140, 200, 300, 600, 1000]
        xedges = [0, 1.4, 2.0, 2.5]
        ybins = np.array(yedges)# + np.diff(yedges)/2
        xbins = np.array(xedges)# + np.diff(xedges)/2
        nameset = ['lf', 'cf', 'bf']
        idx = 0
        for flavor in self.btagSpliter(wp):
            print('drawing 2d hist for',wp ,nameset[idx])
            eff=np.array(unp.nominal_values(flavor[1:,:]))
            #eff=np.flip(eff, axis=0)
            xhist= []
            yhist= []
            weight_2d = []
            for y in range(len(yedges)-1):
                for x in range(len(xedges)-1):
                    xhist.append((xedges[x]+xedges[x+1])/2)
                    yhist.append((yedges[y]+yedges[y+1])/2)
                    weight_2d.append(eff[y][x])
            plt.style.use(hep.style.CMS)
            fig, ax = plt.subplots(1,1,figsize=(10,8))
            hep.cms.label(ax=ax, llabel='Private Work', rlabel=str(year)+' (13 TeV)')
            hist, xbins, ybins, im = plt.hist2d(xhist, yhist, weights=weight_2d, bins=[xbins, ybins], cmap='jet', cmin=0,cmax=1.0)
            for i in range(len(ybins)-1):
                for j in range(len(xbins)-1):
                    ax.text((xbins[j+1] - xbins[j])/2 + xbins[j], (ybins[i+1] - ybins[i])/2 + ybins[i], 
                            np.round(hist.T[i,j],2), 
                            color="w", ha="center", va="top", fontweight="bold")
            # colormap
            plt.colorbar()
            # transpose the array
            ax.set_xlabel(r'Jet $|\eta|$')
            ax.set_ylabel(r'Jet $p_{T}$ [GeV]')
            ax.set_ylim(30,1000)
            ax.set_xlim(0,2.5)
            ax.set_yscale('log')
            plt.savefig(self.outdir+'/btageff_2d_'+wp+'_'+nameset[idx]+'.png')
            plt.close()
            idx += 1

if __name__ == "__main__":
    drawer = Mydrawer(args.i, args.y, args.b)
    #drawer.btagDraw('tight')
    drawer.scaleDraw()
    print("done")
