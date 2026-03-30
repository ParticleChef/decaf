#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from coffea.util import load
import mplhep
import matplotlib.pyplot as plt
import hist

error_opts = {
    'step': 'post',
    'label': 'Stat. ⊕ syst. unc',
    'hatch': '///',
    'facecolor': 'black',
    'alpha': 0.3,
    'edgecolor': (0, 0, 0, 0.3),
    'linewidth': 0
}

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def reduce_to_1d(vals, vars_):
    """
    마지막 축을 x축으로 두고, 나머지 축은 모두 합쳐 1D 배열로 만든다.
    """
    while vals.ndim > 1:
        vals = vals.sum(axis=0)
        if vars_ is not None:
            vars_ = vars_.sum(axis=0)
    return vals, vars_


def get_rebinned_hist(h, key):
    """변수별 rebin 규칙 적용"""
    if key == 'j1pt':
        h = h[{"j1pt": hist.rebin(5)}]
    elif key == 'j2pt':
        h = h[{"j2pt": hist.rebin(5)}]
    elif key == 'fj1pt':
        h = h[{"fj1pt": hist.rebin(5)}]
    elif key == 'fj1TvsQCD':
        h = h[{"fj1TvsQCD": hist.rebin(10)}]
    elif key == 'metpt_10GeVbins':
        h = h[{"metpt_10GeVbins": hist.rebin(5)}]
    elif key == 'recoilpt':
        h = h[{"recoilpt": hist.rebin(5)}]
    return h


def get_hist_safe(container, key, process, region, systematic):
    try:
        h = container[key][process][{'region': region, 'systematic': systematic}]
        h = get_rebinned_hist(h, key)
        return h
    except Exception:
        return None


def get_data_process(region):
    if 'GCR' in region or 'DY2E' in region:
        return 'EGamma'
    elif 'DY2M' in region:
        return 'Muon'
    else:
        return 'JetMET'


def get_pretty_label(mc):
    if mc == 'Z (inv)':
        return r'Z$\rightarrow$$\nu$$\bar{\nu}$'
    elif mc == 'W (lnu)':
        return r'W$\rightarrow$$\ell$$\nu$'
    elif mc == 'TT':
        return r't$\bar{t}$'
    else:
        return mc


# ----------------------------------------------------------------------
# Load histograms
# ----------------------------------------------------------------------
myhist = load('hists/stop_new2024.scaled')
bkg = myhist['bkg']
data = myhist['data']
sig = myhist['sig']

# ----------------------------------------------------------------------
# Plot config
# ----------------------------------------------------------------------
nominal_syst = 'nominal'

# nuisance별 Up/Down systematic 이름
syst_pairs = {
    'pileup': ('pileupUp', 'pileupDown'),
    'electron_id': ('electron_idUp', 'electron_idDown'),
    'muon_id': ('muon_idUp', 'muon_idDown'),
}

lumi_unc = 0.016  # 1.6%

stacks = ['VV', 'Single Top', 'TT', 'DY', 'Gamma + Jets', 'W (lnu)', 'Z (inv)', 'QCD Multijet']
colors = ['#6b705c', '#8e7dbe', '#99c1b9', 'cyan', 'purple', '#f1e3d3', '#f2d0a9', '#d88c9a']
signals = ['SMS-2Stop-Par-mStop-1000', 'SMS-2Stop-Par-mStop-1500', 'SMS-2Stop-Par-mStop-600']

signal_labels = [
    r'$m_{\tilde{t}}$ = 1000 GeV',
    r'$m_{\tilde{t}}$ = 1500 GeV',
    r'$m_{\tilde{t}}$ = 600 GeV'
]
signal_ltypes = ['solid', 'dashed', 'dotted']
signal_colors = ['red', 'darkred', 'maroon']

regions_to_plot = [
    'cat1_preselection',
    'cat2_LLCR_highDeltaM',
    'cat3_QCDCR_highDeltaM',
    'cat4_GCR_highDeltaM',
    'cat5_DY2E_highDeltaM',
    'cat6_DY2M_highDeltaM',
    'cat7_SR_highDeltaM'
]

os.makedirs('plots', exist_ok=True)

print(bkg['metpt']['W (lnu)'])

# ----------------------------------------------------------------------
# Print yields for checking
# ----------------------------------------------------------------------
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

print("Systematics keys in bkg['metpt']['W (lnu)']:", bkg['metpt']['W (lnu)'])

# ----------------------------------------------------------------------
# Main plotting loop
# ----------------------------------------------------------------------
for region in regions_to_plot:
    print(f"Plotting for region: {region}")

    for key in bkg.keys():
        if 'sumw' in key:
            continue
        elif 'template' in key:
            continue
        elif 'nPV' in key:
            continue

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
        # 1) MC hist collection
        # --------------------------------------------------------------
        mc_vals_list = []
        mc_vars_list = []
        mc_labels = []
        bins = None

        mc_vals_syst_up = {name: [] for name in syst_pairs}
        mc_vals_syst_down = {name: [] for name in syst_pairs}

        for mc in stacks:
            if mc not in bkg[key]:
                continue

            h_nom = get_hist_safe(bkg, key, mc, region, nominal_syst)
            if h_nom is None:
                continue

            if bins is None:
                bins = h_nom.axes[-1].edges

            vals_nom = h_nom.values()
            vars_nom = h_nom.variances()

            vals_nom_1d, vars_nom_1d = reduce_to_1d(vals_nom, vars_nom)
            mc_vals_list.append(vals_nom_1d)
            mc_vars_list.append(vars_nom_1d)
            mc_labels.append(get_pretty_label(mc))

            for syst_name, (syst_up_name, syst_down_name) in syst_pairs.items():
                h_up = get_hist_safe(bkg, key, mc, region, syst_up_name)
                if h_up is not None:
                    vals_up_1d, _ = reduce_to_1d(h_up.values(), None)
                else:
                    vals_up_1d = vals_nom_1d.copy()

                h_down = get_hist_safe(bkg, key, mc, region, syst_down_name)
                if h_down is not None:
                    vals_down_1d, _ = reduce_to_1d(h_down.values(), None)
                else:
                    vals_down_1d = vals_nom_1d.copy()

                mc_vals_syst_up[syst_name].append(vals_up_1d)
                mc_vals_syst_down[syst_name].append(vals_down_1d)

        if len(mc_vals_list) == 0 or bins is None:
            plt.close()
            continue

        mc_vals_arr = np.array(mc_vals_list)
        mc_vars_arr = np.array(mc_vars_list)

        # --------------------------------------------------------------
        # 2) Total background + uncertainties
        # --------------------------------------------------------------
        total_bkg = mc_vals_arr.sum(axis=0)
        total_bkg_var = mc_vars_arr.sum(axis=0)
        total_bkg_err = np.sqrt(total_bkg_var)

        # nuisance별 envelope -> nuisance끼리 quadrature
        sys_up_sq = np.zeros_like(total_bkg, dtype=float)
        sys_down_sq = np.zeros_like(total_bkg, dtype=float)

        total_bkg_syst_up = {}
        total_bkg_syst_down = {}
        per_syst_bin_up = {}
        per_syst_bin_down = {}

        for syst_name in syst_pairs:
            mc_vals_up_arr = np.array(mc_vals_syst_up[syst_name])
            mc_vals_down_arr = np.array(mc_vals_syst_down[syst_name])

            total_up = mc_vals_up_arr.sum(axis=0)
            total_down = mc_vals_down_arr.sum(axis=0)

            total_bkg_syst_up[syst_name] = total_up
            total_bkg_syst_down[syst_name] = total_down

            delta_up_var = total_up - total_bkg
            delta_down_var = total_down - total_bkg

            this_sys_up = np.maximum.reduce([
                delta_up_var,
                delta_down_var,
                np.zeros_like(total_bkg)
            ])
            this_sys_down = np.maximum.reduce([
                -delta_up_var,
                -delta_down_var,
                np.zeros_like(total_bkg)
            ])

            per_syst_bin_up[syst_name] = this_sys_up
            per_syst_bin_down[syst_name] = this_sys_down

            sys_up_sq += this_sys_up**2
            sys_down_sq += this_sys_down**2

        sys_up = np.sqrt(sys_up_sq)
        sys_down = np.sqrt(sys_down_sq)

        # lumi: normalization only
        lumi_err = total_bkg * lumi_unc

        # final total band
        tot_err_up = np.sqrt(total_bkg_err**2 + sys_up**2 + lumi_err**2)
        tot_err_down = np.sqrt(total_bkg_err**2 + sys_down**2 + lumi_err**2)

        unc_low = total_bkg - tot_err_down
        unc_up = total_bkg + tot_err_up
        unc_low = np.maximum(unc_low, 1e-10)

        total_bkg_forratio = np.append(total_bkg, total_bkg[-1])
        unc_low_ratio = np.append(unc_low, unc_low[-1])
        unc_up_ratio = np.append(unc_up, unc_up[-1])

        # --------------------------------------------------------------
        # Debug print (metpt only)
        # --------------------------------------------------------------
        if key == 'metpt':
            print(f"\n[DEBUG] region={region}, variable={key}")

            nom_yield = total_bkg.sum()
            print(f"  Total yield nominal = {nom_yield:.6f}")

            for syst_name in syst_pairs:
                up_yield = total_bkg_syst_up[syst_name].sum()
                down_yield = total_bkg_syst_down[syst_name].sum()

                rel_up_yield = 100.0 * (up_yield - nom_yield) / nom_yield if nom_yield != 0 else 0.0
                rel_down_yield = 100.0 * (down_yield - nom_yield) / nom_yield if nom_yield != 0 else 0.0

                print(f"  {syst_name:12s} Up   = {up_yield:.6f}  ({rel_up_yield:+.3f}%)")
                print(f"  {syst_name:12s} Down = {down_yield:.6f}  ({rel_down_yield:+.3f}%)")

                with np.errstate(divide='ignore', invalid='ignore'):
                    rel_diff_up = np.divide(
                        total_bkg_syst_up[syst_name] - total_bkg,
                        total_bkg,
                        out=np.zeros_like(total_bkg, dtype=float),
                        where=(total_bkg != 0)
                    ) * 100.0

                    rel_diff_down = np.divide(
                        total_bkg_syst_down[syst_name] - total_bkg,
                        total_bkg,
                        out=np.zeros_like(total_bkg, dtype=float),
                        where=(total_bkg != 0)
                    ) * 100.0

                max_abs_up = np.max(np.abs(total_bkg_syst_up[syst_name] - total_bkg))
                max_abs_down = np.max(np.abs(total_bkg_syst_down[syst_name] - total_bkg))
                max_rel_up = np.max(np.abs(rel_diff_up))
                max_rel_down = np.max(np.abs(rel_diff_down))

                print(f"    Max |Up - Nominal| per bin   = {max_abs_up:.6f}")
                print(f"    Max |Down - Nominal| per bin = {max_abs_down:.6f}")
                print(f"    Max |Up - Nominal| / Nominal per bin   = {max_rel_up:.3f}%")
                print(f"    Max |Down - Nominal| / Nominal per bin = {max_rel_down:.3f}%")

            print(f"  Lumi uncertainty = {100.0 * lumi_unc:.3f}%")

            # total yield uncertainty
            yield_sys_up_sq = 0.0
            yield_sys_down_sq = 0.0

            for syst_name in syst_pairs:
                up_yield = total_bkg_syst_up[syst_name].sum()
                down_yield = total_bkg_syst_down[syst_name].sum()

                delta_up = up_yield - nom_yield
                delta_down = down_yield - nom_yield

                this_up = max(delta_up, delta_down, 0.0)
                this_down = max(-delta_up, -delta_down, 0.0)

                yield_sys_up_sq += this_up**2
                yield_sys_down_sq += this_down**2

            stat_yield_err = np.sqrt(np.sum(total_bkg_var))
            lumi_yield_err = nom_yield * lumi_unc

            total_yield_err_up = np.sqrt(stat_yield_err**2 + yield_sys_up_sq + lumi_yield_err**2)
            total_yield_err_down = np.sqrt(stat_yield_err**2 + yield_sys_down_sq + lumi_yield_err**2)

            rel_total_yield_err_up = 100.0 * total_yield_err_up / nom_yield if nom_yield != 0 else 0.0
            rel_total_yield_err_down = 100.0 * total_yield_err_down / nom_yield if nom_yield != 0 else 0.0

            print(
                f"  Total yield uncertainty       = "
                f"+{total_yield_err_up:.6f} ({rel_total_yield_err_up:.3f}%), "
                f"-{total_yield_err_down:.6f} ({rel_total_yield_err_down:.3f}%)"
            )

            # band-like combined uncertainty
            total_err_up_yield = np.sqrt(np.sum(tot_err_up**2))
            total_err_down_yield = np.sqrt(np.sum(tot_err_down**2))

            rel_total_err_up_yield = 100.0 * total_err_up_yield / nom_yield if nom_yield != 0 else 0.0
            rel_total_err_down_yield = 100.0 * total_err_down_yield / nom_yield if nom_yield != 0 else 0.0

            print(
                f"  Band-like combined uncertainty= "
                f"+{total_err_up_yield:.6f} ({rel_total_err_up_yield:.3f}%), "
                f"-{total_err_down_yield:.6f} ({rel_total_err_down_yield:.3f}%)"
            )

        # --------------------------------------------------------------
        # 3) Stacked MC
        # --------------------------------------------------------------
        bottom = np.zeros_like(bins[:-1], dtype=float)
        stack_tops = []

        for vals_1d, lab, col in zip(mc_vals_list, mc_labels, colors):
            top = bottom + vals_1d

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

        for top in stack_tops:
            y = np.append(top, top[-1])
            ax.step(
                bins,
                y,
                where='post',
                color='black',
                linewidth=0.8,
            )

        # --------------------------------------------------------------
        # 4) Total uncertainty band
        # --------------------------------------------------------------
        ax.fill_between(bins, unc_low_ratio, unc_up_ratio, **error_opts)

        # --------------------------------------------------------------
        # 5) Signal draw
        # --------------------------------------------------------------
        if 'CR' not in region and 'DY' not in region:
            for i, signal in enumerate(signals):
                if signal not in bkg[key]:
                    continue

                h_sig = get_hist_safe(bkg, key, signal, region, nominal_syst)
                if h_sig is None:
                    continue

                    vals_sig = h_sig.values()
                vals_sig = h_sig.values()
                vals_sig_1d, _ = reduce_to_1d(vals_sig, None)
                sig_bins = h_sig.axes[-1].edges

                ax.step(
                    sig_bins,
                    np.append(vals_sig_1d, vals_sig_1d[-1]),
                    where='post',
                    label=signal_labels[i],
                    color=signal_colors[i],
                    linestyle=signal_ltypes[i],
                    linewidth=2.5,
                )

        # --------------------------------------------------------------
        # 6) Data
        # --------------------------------------------------------------
        data_process = get_data_process(region)
        try:
            h_data = data[key][data_process][{'region': region, 'systematic': nominal_syst}]
            h_data = get_rebinned_hist(h_data, key)
        except Exception:
            plt.close()
            continue

        bins_data = h_data.axes[-1].edges
        centers = h_data.axes[-1].centers

        vals_data = h_data.values()
        vals_data_1d, _ = reduce_to_1d(vals_data, None)

        bins = bins_data
        yerr_data = np.sqrt(vals_data_1d)

        if 'SR' not in region:
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

        if 'CR' in region or 'DY' in region:
            ax.set_ylim(0.1, 1e+8)
        elif 'SR' in region:
            ax.set_ylim(0.1, 1e+8)
        else:
            ax.set_ylim(1, 1e+10)

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
        if 'mll' in key:
            ax.set_xlim(50, 250)
        if 'pll' in key:
            ax.set_xlim(200, 1000)

        ax.grid(True, which='both', axis='y', ls='--', lw=0.5)

        # --------------------------------------------------------------
        # 7) Ratio plot
        # --------------------------------------------------------------
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.divide(
                vals_data_1d,
                total_bkg,
                out=np.zeros_like(vals_data_1d, dtype=float),
                where=(total_bkg != 0)
            )
            ratio_err = np.divide(
                yerr_data,
                total_bkg,
                out=np.zeros_like(yerr_data, dtype=float),
                where=(total_bkg != 0)
            )
            ratio_err = np.where(ratio_err < 0, 0, ratio_err)

        if 'SR' not in region:
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
                unc_low_ratio,
                total_bkg_forratio,
                out=np.ones_like(unc_low_ratio, dtype=float),
                where=(total_bkg_forratio != 0)
            )
            band_up = np.divide(
                unc_up_ratio,
                total_bkg_forratio,
                out=np.ones_like(unc_up_ratio, dtype=float),
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
            rax.set_xlabel(r"$E\!\!\!/_{T}\ \mathrm{(GeV)}$")
        elif 'recoilpt' in key:
            rax.set_xlabel(r"$U\!\!\!/_{T}\ \mathrm{(GeV)}$")
        elif 'ht' in key and 'tight' not in key:
            rax.set_xlabel('H$_{T}$ (GeV)')
            ax.set_xlim(300, 1500)
            rax.set_xlim(300, 1500)
        elif 'nJet' in key:
            rax.set_xlabel('Number of Jets')
            ax.set_xlim(2, 10)
            rax.set_xlim(2, 10)
        elif 'nElectron' in key:
            rax.set_xlabel('Number of Electrons')
        elif 'nMuon' in key:
            rax.set_xlabel('Number of Muons')
        elif 'nb' in key:
            rax.set_xlabel('Number of b jets')
            ax.set_xlim(0, 5)
            rax.set_xlim(0, 5)
        elif key == 'j1pt':
            rax.set_xlabel('Leading jet $p_{T}$ [GeV]')
            ax.set_xlim(100, 600)
            rax.set_xlim(100, 600)
        elif key == 'j2pt':
            rax.set_xlabel('Subleading jet $p_{T}$ [GeV]')
            ax.set_xlim(50, 600)
            rax.set_xlim(50, 600)

        xmin, xmax = ax.get_xlim()
        rax.set_xlim(xmin, xmax)

        rax.set_ylabel('Data/MC')
        rax.set_ylim(0, 2)
        rax.set_yticks([0, 0.5, 1, 1.5, 2])
        rax.grid(True, which='both', axis='y', ls='--', lw=0.5)

        plt.savefig(f'plots/{key}_{region}.png')
        plt.close()

print("Done.")