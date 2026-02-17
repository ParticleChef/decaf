#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from coffea.util import load
import mplhep
import matplotlib.pyplot as plt
from scipy.special import betaincinv  # type: ignore
import hist

error_opts = {
    'step': 'post',
    'label': 'Stat. unc',
    'hatch': '///',
    'facecolor': 'black',
    'alpha': 0.3,
    'edgecolor': (0, 0, 0, 0.3),
    'linewidth': 0
}

# ----------------------------------------------------------------------
# Load histograms
# ----------------------------------------------------------------------
myhist = load('hists/stop_new2024.scaled')
bkg = myhist['bkg']
data = myhist['data']
sig = myhist['sig']
## signal
#print(sig)
#print(data)
#print(sig['metpt'].keys())


# ----------------------------------------------------------------------
# Plot config
# ----------------------------------------------------------------------

systematics = 'nominal'

# MC stack 순서/색
#stacks = ['W (lnu)', 'QCD Multijet', 'TT', 'Single Top']
#colors = ['green', 'lightblue', 'yellow', 'orange']
#stacks = ['VV','Single Top', 'TT',r"$\gamma$ + Jets", 'W (lnu)','Z (inv)', 'QCD Multijet']
stacks = ['VV','Single Top', 'TT','DY','Gamma + Jets', 'W (lnu)','Z (inv)', 'QCD Multijet']
#colors = ['','orange', 'yellow', 'green','salmon', 'lightblue']
colors = ['#6b705c','#8e7dbe', '#99c1b9','cyan','purple' ,'#f1e3d3', '#f2d0a9', '#d88c9a']
#stacks = ['Single Top', 'TT', 'QCD Multijet']
#colors = ['orange', 'yellow', 'lightblue']
signals = ['SMS-2Stop-Par-mStop-1000', 'SMS-2Stop-Par-mStop-1500', 'SMS-2Stop-Par-mStop-600']

print(bkg['metpt']['W (lnu)'])

### print yields for checking
total_bkg_yield = 0
for process in stacks:
    if process in bkg['nMuon']:
        h = bkg['nMuon'][process][{
            'region': 'cat1_preselection',
            'systematic': 'nominal'
        }]
        yield_ = h.values().sum()
        total_bkg_yield += yield_
        print(f'Yield for {process}: {yield_}')
print(f'Total background yield: {total_bkg_yield}')
### print data yield
h_data = data['nMuon']['JetMET'][{
    'region': 'cat1_preselection',
    'systematic': 'nominal'
}]
data_yield = h_data.values().sum()
print(f'Data yield: {data_yield}')

for signal in signals:
    if signal in bkg['nMuon']:
        h_sig = bkg['nMuon'][signal][{
            'region': 'cat1_preselection',
            'systematic': 'nominal'
        }]
        sig_yield = h_sig.values().sum()
        print(f'Yield for {signal}: {sig_yield}')

def reduce_to_1d(vals, vars_):
    """
    vals, vars_: 같은 shape의 ndarray.
    마지막 축을 x축(플로팅용)으로 두고, 나머지 축들은 모두 합쳐서 1D로 만든다.
    """
    while vals.ndim > 1:
        vals = vals.sum(axis=0)
        if vars_ is not None:
            vars_ = vars_.sum(axis=0)
    return vals, vars_


# ----------------------------------------------------------------------
# Loop over variables in bkg (e.g. metpt, metphi, nJet, ...)
# ----------------------------------------------------------------------
for regions in ['cat6_DY2M_highDeltaM']:#['cat1_preselection', 'cat2_LLCR_highDeltaM', 'cat3_QCDCR_highDeltaM', 'cat4_GCR_highDeltaM', 'cat5_DY2E_highDeltaM', 'cat6_DY2M_highDeltaM', 'cat7_SR_highDeltaM']:
    print(f"Plotting for region: {regions}")
    for key in bkg.keys():
        if 'sumw' in key:
            continue
        elif 'template' in key:
            continue
        elif 'nPV' in key:
            continue
        else:
            #print("Plotting ", key)

            plt.style.use(mplhep.style.CMS)
            fig, (ax, rax) = plt.subplots(
                nrows=2,
                ncols=1,
                figsize=(9, 9),
                gridspec_kw={"height_ratios": (3, 1)},
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.07)
            mplhep.cms.label(
                ax=ax,
                llabel='Work in progress',
                rlabel='108.95 fb$^{-1}$ (13.6 TeV)',
                fontsize=24
            )

            # --------------------------------------------------------------
            # 1) MC 히스토그램 수집 (numpy 배열만 사용)
            # --------------------------------------------------------------
            mc_vals_list = []   # 각 MC의 1D 값
            mc_vars_list = []   # 각 MC의 1D 분산
            mc_labels = []
            bins = None

            for mc in stacks:
                if mc not in bkg[key]:
                    continue
                try:
                    h = bkg[key][mc][{'region': regions, 'systematic': systematics}]
                except:
                    continue
                if key == 'j1pt':
                    # rebin for j1pt
                    h = h[{"j1pt": hist.rebin(5)}]
                elif key == 'j2pt':
                    # rebin for j2pt
                    h = h[{"j2pt": hist.rebin(5)}]
                elif key == 'fj1pt':
                    # rebin for fj1pt
                    h = h[{"fj1pt": hist.rebin(5)}]
                elif key == 'fj1TvsQCD':
                    # rebin for fj1TvsQCD
                    h = h[{"fj1TvsQCD": hist.rebin(10)}]
                elif key == 'metpt_10GeVbins':
                    # rebin for metpt_10GeVbins
                    h = h[{"metpt_10GeVbins": hist.rebin(5)}]
                elif key == 'recoilpt':
                    # rebin for recoilpt
                    h = h[{"recoilpt": hist.rebin(5)}]

                # 모든 축 중 마지막 축만 x축으로 두고, 나머지는 reduce_to_1d에서 합침
                if bins is None:
                    bins = h.axes[-1].edges
                else:
                    # bin 같다고 가정 (다르면 여기서 체크 가능)
                    pass

                vals = h.values()
                vars_ = h.variances()

                vals_1d, vars_1d = reduce_to_1d(vals, vars_)
                mc_vals_list.append(vals_1d)
                mc_vars_list.append(vars_1d)

                # legend label
                if 'Z (inv)' in mc:
                    lab = r'Z$\rightarrow$$\nu$$\bar{\nu}$'
                elif 'W (lnu)' in mc:
                    lab = r'W$\rightarrow$$\ell$$\nu$'
                elif 'TT' in mc:
                    lab = r't$\bar{t}$'
                else:
                    lab = mc
                mc_labels.append(lab)

            # MC가 하나도 없으면 스킵
            if len(mc_vals_list) == 0 or bins is None:
                plt.close()
                continue

            mc_vals_arr = np.array(mc_vals_list)       # shape: (nMC, nbins)
            mc_vars_arr = np.array(mc_vars_list)       # shape: (nMC, nbins)

            # --------------------------------------------------------------
            # 2) total background + stat. uncertainty 계산 (배열 기반)
            # --------------------------------------------------------------
            total_bkg = mc_vals_arr.sum(axis=0)                # nbins
            total_bkg_var = mc_vars_arr.sum(axis=0)            # nbins
            total_bkg_err = np.sqrt(total_bkg_var)

            unc_low = total_bkg - total_bkg_err
            unc_up = total_bkg + total_bkg_err

            # ratio plot용으로 맨 마지막 bin 하나 더 붙임
            total_bkg_forratio = np.append(total_bkg, total_bkg[-1])
            unc_low_ratio = np.append(unc_low, unc_low[-1])
            unc_up_ratio = np.append(unc_up,   unc_up[-1])

            # --------------------------------------------------------------
            # 3) Stacked MC 그리기 (bottom 누적 방식)
            # --------------------------------------------------------------
            bottom = np.zeros_like(bins[:-1])

            stack_tops = []  # 각 층의 윗부분 값 (경계선용)

            for vals_1d, lab, col in zip(mc_vals_list, mc_labels, colors):
                top = bottom + vals_1d  # 이 프로세스까지 쌓였을 때의 높이

                # 먼저 색 채우기 (테두리는 없음)
                ax.bar(
                    bins[:-1],
                    vals_1d,
                    width=np.diff(bins),
                    bottom=bottom,
                    align='edge',
                    label=lab,
                    color=col,
                    edgecolor='none',
                    linewidth=0,
                )

                stack_tops.append(top.copy())
                bottom = top

            # 각 층의 윗면에 얇은 검은 실선 그리기
            for top in stack_tops:
                y = np.append(top, top[-1])  # step='post' 위해 마지막 bin 하나 복제
                ax.step(
                    bins,
                    y,
                    where='post',
                    color='black',
                    linewidth=0.8,
                )

            # --------------------------------------------------------------
            # 4) MC stat. uncertainty band (total bkg)
            #    → x: bins[:-1] (24), y: unc_low/up (24)로 길이 맞추기
            # --------------------------------------------------------------
            ax.fill_between(bins, unc_low_ratio, unc_up_ratio, **error_opts)

            labs = [r'$m_{\tilde{t}}$ = 1000 GeV', r'$m_{\tilde{t}}$ = 1500 GeV', r'$m_{\tilde{t}}$ = 600 GeV']
            ltypes = ['solid', 'dashed', 'dotted']
            scolors = ['red', 'darkred', 'maroon']

            ## signal draw
            if 'CR' in regions or 'DY' in regions:
                pass
            else:
                for signal in signals:
                    if signal not in bkg[key]:
                        continue
                    try:
                        h_sig = bkg[key][signal][{'region': regions, 'systematic': systematics}]
                    except:
                        continue
                    if key == 'j1pt':
                        # rebin for j1pt
                        h_sig = h_sig[{"j1pt": hist.rebin(5)}]
                    elif key == 'j2pt':
                        # rebin for j2pt
                        h_sig = h_sig[{"j2pt": hist.rebin(5)}]
                    elif key == 'fj1pt':
                        # rebin for fj1pt
                        h_sig = h_sig[{"fj1pt": hist.rebin(5)}]
                    elif key == 'fj1TvsQCD':
                        # rebin for fj1TvsQCD
                        h_sig = h_sig[{"fj1TvsQCD": hist.rebin(10)}]
                    elif key == 'metpt_10GeVbins':
                        # rebin for metpt_10GeVbins
                        h_sig = h_sig[{"metpt_10GeVbins": hist.rebin(5)}]
                    elif key == 'recoilpt':
                        # rebin for recoilpt
                        h_sig = h_sig[{"recoilpt": hist.rebin(5)}]


                    vals_sig = h_sig.values()
                    ## scaling to original cross-section
                    #if '600' in signal:
                    #    vals_sig *= (2.560e-01 / 10)
                    #elif '1000' in signal:
                    #    vals_sig *= (9.123e-03 / 10)
                    #elif '1500' in signal:
                    #    vals_sig *= (3.912e-04 / 10)
                    vals_sig_1d, _ = reduce_to_1d(vals_sig, None)

                    # bins는 MC 기준으로 통일 (data와 같다고 가정)
                    bins = h_sig.axes[-1].edges
                    centers = h_sig.axes[-1].centers

                    ax.step(
                        bins,
                        np.append(vals_sig_1d, vals_sig_1d[-1]),
                        where='post',
                        label=labs[signals.index(signal)],
                        color=scolors[signals.index(signal)],
                        linestyle=ltypes[signals.index(signal)],
                        linewidth=2.5,
                    )


            # --------------------------------------------------------------
            # 5) Data points + error bars (배열로 reduce)
            # --------------------------------------------------------------
            if 'GCR' in regions or 'DY2E' in regions:
                h_data = data[key]['EGamma'][{'region': regions, 'systematic': systematics}]
            elif 'DY2M' in regions:
                h_data = data[key]['Muon'][{'region': regions, 'systematic': systematics}]
            else:
                h_data = data[key]['JetMET'][{'region': regions, 'systematic': systematics}]
            if key == 'j1pt':
                # rebin for j1pt
                h_data = h_data[{"j1pt": hist.rebin(5)}]
            elif key == 'j2pt':
                # rebin for j2pt
                h_data = h_data[{"j2pt": hist.rebin(5)}]
            elif key == 'fj1pt':
                # rebin for fj1pt
                h_data = h_data[{"fj1pt": hist.rebin(5)}]
            elif key == 'fj1TvsQCD':
                # rebin for fj1TvsQCD
                h_data = h_data[{"fj1TvsQCD": hist.rebin(10)}]
            elif key == 'metpt_10GeVbins':
                # rebin for metpt_10GeVbins
                h_data = h_data[{"metpt_10GeVbins": hist.rebin(5)}]
            elif key == 'recoilpt':
                # rebin for recoilpt
                h_data = h_data[{"recoilpt": hist.rebin(5)}]


            bins_data = h_data.axes[-1].edges
            centers = h_data.axes[-1].centers

            vals_data = h_data.values()
            vals_data_1d, _ = reduce_to_1d(vals_data, None)

            # bins는 data 기준으로 통일 (MC와 같다고 가정)
            bins = bins_data
            centers = h_data.axes[-1].centers

            yerr_data = np.sqrt(vals_data_1d)

            if 'SR' in regions:
                pass
            else:
                ax.errorbar(
                    centers,
                    vals_data_1d,
                    xerr=np.diff(bins) / 2,
                    yerr=yerr_data,
                    fmt='o',
                    label='Data',
                    markersize=8,
                    capsize=5,
                    capthick=1,
                    color='black',
                )

            ax.set_ylabel('Events')
            ax.legend(ncol=3, fontsize=14, loc='upper right')
            ax.set_yscale('log')
            ax.set_yticks([1, 10, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9, 1e+10])
            if 'CR' in regions or 'DY' in regions:
                ax.set_ylim(0.1, 1e+8)
            elif 'SR' in regions:
                ax.set_ylim(0.1, 1e+8)
            else:
                ax.set_ylim(1, 1e+10)

            # 변수에 따라 x-range 조정 (원래 코드 로직 유지)
            if 'metpt' in key:
                ax.set_xlim(250, 800)
            if 'recoilpt' in key:
                ax.set_xlim(200, 800)
            if 'eta' in key:
                ax.set_xlim(-3.0, 3.0)
            if 'phi' in key:
                ax.set_xlim(-3.2, 3.2)
            if 'nElectron' in key or 'nMuon' in key:
                ax.set_xlim(0, 4)
            if 'nJet' in key:
                ax.set_xlim(0, 10)
            if 'nfj' in key:
                ax.set_xlim(0, 6)
            if 'vs' in key:
                ax.set_xlim(0, 1)
            if key == 'fj1pt':
                ax.set_xlim(200, 1000)

            ax.grid(True, which='both', axis='y', ls='--', lw=0.5)

            # --------------------------------------------------------------
            # 6) Ratio plot (Data/MC) + MC uncertainty band
            # --------------------------------------------------------------
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.divide(
                    vals_data_1d, total_bkg,
                    out=np.zeros_like(vals_data_1d),
                    where=(total_bkg != 0)
                )
                ratio_err = np.divide(
                    yerr_data, total_bkg,
                    out=np.zeros_like(yerr_data),
                    where=(total_bkg != 0)
                )
                ## Negative value to zero
                ratio_err = np.where(ratio_err < 0, 0, ratio_err)
            if 'SR' in regions:
                pass
            else:
                rax.errorbar(
                    centers,
                    ratio,
                    xerr=np.diff(bins) / 2,
                    yerr=ratio_err,
                    fmt='o',
                    color='black',
                    markersize=8,
                    capsize=5,
                    capthick=1,
                )

            with np.errstate(divide='ignore', invalid='ignore'):
                band_low = np.divide(
                    unc_low_ratio, total_bkg_forratio,
                    out=np.ones_like(unc_low_ratio),
                    where=(total_bkg_forratio != 0)
                )
                band_up = np.divide(
                    unc_up_ratio, total_bkg_forratio,
                    out=np.ones_like(unc_up_ratio),
                    where=(total_bkg_forratio != 0)
                )

            rax.fill_between(
                bins,
                band_low,
                band_up,
                **error_opts
            )

            rax.axhline(1.0, color='black', linestyle='--', linewidth=1)
            rax.set_xlabel(key)

            if 'metpt' in key:
                #rax.set_xlabel('Missing $p_{T}$ [GeV]')
                rax.set_xlabel(r"$E\!\!\!/_{T}\ \mathrm{(GeV)}$")
            elif 'recoilpt' in key:
                rax.set_xlabel(r"$U\!\!\!/_{T}\ \mathrm{(GeV)}$")
            elif 'ht' in key and not 'tight' in key:
                rax.set_xlabel('H$_{T}$ (GeV)')
                ax.set_xlim(300, 1500)
                rax.set_xlim(300, 1500)
            elif 'nJet' in key:
                rax.set_xlabel('Number of Jets')
                ax.set_xlim(2, 10)
            elif 'nElectron' in key:
                rax.set_xlabel('Number of Electrons')
            elif 'nMuon' in key:
                rax.set_xlabel('Number of Muons')
            elif 'nb' in key:
                rax.set_xlabel('Number of b jets')
                ax.set_xlim(0, 5)
            elif key == 'j1pt':
                rax.set_xlabel('Leading jet $p_{T}$ [GeV]')
                ax.set_xlim(100, 600)
            elif key == 'j2pt':
                rax.set_xlabel('Subleading jet $p_{T}$ [GeV]')
                ax.set_xlim(50, 600)
                
            
            rax.set_ylabel('Data/MC')
            rax.set_ylim(0, 2)
            rax.set_yticks([0, 0.5, 1, 1.5, 2])
            rax.grid(True, which='both', axis='y', ls='--', lw=0.5)

            plt.savefig(f'plots/{key}_{regions}.png')
            plt.close()