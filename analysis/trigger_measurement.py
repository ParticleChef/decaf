import uproot
import numpy as np
import awkward as ak
from coffea.util import load, save
import mplhep
import matplotlib.pyplot as plt
import hist
from scipy.special import betaincinv # type: ignore

def ClopperPearson(total, passed, level=0.68):
    """
    Clopper-Pearson method for calculating confidence intervals.
    :param total: Array of total events
    :param passed: Array of events that passed the selection
    :param level: Confidence level (0.68 for 1 sigma, 0.95 for 2 sigma, etc.)
    :return: Lower and upper bounds of the confidence interval
    """
    alpha = 1 - level
    lower_bound = betaincinv(passed, total - passed + 1, alpha / 2)
    upper_bound = betaincinv(passed + 1, total - passed, 1 - alpha / 2)
    # nan to 0 for lower bound and 1 for upper bound
    lower_bound = np.where(np.isnan(lower_bound), 0, lower_bound)
    upper_bound = np.where(np.isnan(upper_bound), 1, upper_bound)
    return lower_bound, upper_bound

# Load the merged or scaled histograms
def get_trigger_efficiency(scaled):
    myhist = load(scaled)
    h = myhist['data']['metpt']['EGamma'][{'region': 'cat1_preselection', 'systematic': 'nominal'}]
    h_tot = (h[{'signal_trigger': 0}] + h[{'signal_trigger': 1}])
    h_pass = (h[{'signal_trigger': 1}])
    bins = h_tot.axes[0].edges
    centers  = h_tot.axes[0].centers
    eff = np.divide(h_pass.values(), h_tot.values(), out=np.zeros_like(h_pass.values()), where=h_tot.values() != 0)
    print("Efficiency:", eff)
    err_low, err_up = ClopperPearson(h_tot.values(), h_pass.values(), level=0.68)
    return bins, centers, eff, err_low, err_up


EGamma_2022 = get_trigger_efficiency('hists/egamma.scaled')
EGamma_2023 = get_trigger_efficiency('hists/2023_EGamma.scaled')
EGamma_2024 = get_trigger_efficiency('hists/stop_2024.scaled')

# Plotting the efficiency with error bars
plt.figure(figsize=(8, 3))
plt.style.use(mplhep.style.CMS)
mplhep.cms.label(llabel='Work in progress', rlabel='(13.6 TeV)')
'''
plt.errorbar(
    centers, eff,
    xerr = np.diff(bins) / 2,
    yerr=[np.abs(eff - err_low), np.abs(eff - err_up)],
    fmt='o', label='2022 EGamma (C,D,E,F,G)',
    markersize=8, capsize=5, capthick=1,
    color='black',
)
'''
# 2022
plt.errorbar(
    EGamma_2022[1], EGamma_2022[2],
    xerr = np.diff(EGamma_2022[0]) / 2,
    yerr=[np.abs(EGamma_2022[2] - EGamma_2022[3]), np.abs(EGamma_2022[4] - EGamma_2022[2])],
    fmt='o', label='2022 EGamma (C,D,E,F,G)\n 34.65 fb$^{-1}$',
    markersize=8, capsize=5, capthick=1,
    color='black',
)
# 2023
plt.errorbar(
    EGamma_2023[1], EGamma_2023[2],
    xerr = np.diff(EGamma_2023[0]) / 2,
    yerr=[np.abs(EGamma_2023[2] - EGamma_2023[3]), np.abs(EGamma_2023[4] - EGamma_2023[2])],
    fmt='s', label='2023 EGamma (C,D)\n 27.25 fb$^{-1}$',
    markersize=8, capsize=5, capthick=1,
    color='red',
)
# 2024
plt.errorbar(
    EGamma_2024[1], EGamma_2024[2],
    xerr = np.diff(EGamma_2024[0]) / 2,
    yerr=[np.abs(EGamma_2024[2] - EGamma_2024[3]), np.abs(EGamma_2024[4] - EGamma_2024[2])],
    fmt='^', label='2024 EGamma (C,D,E,F,G,H,I)\n 108.95 fb$^{-1}$',
    markersize=8, capsize=5, capthick=1,
    color='blue',
)

plt.xlabel('$E^{miss}_{T}$ (GeV)')
plt.ylabel('Trigger Efficiency')

plt.xticks(np.arange(0, 900, 100))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.xlim(100,800)
plt.ylim(0, 1.01)
plt.grid()
plt.legend(loc='lower right', fontsize=20)
plt.tight_layout()
plt.savefig('trigger_efficiency2.png')