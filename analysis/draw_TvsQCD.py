import numpy as np
from coffea.util import load
import mplhep
import matplotlib.pyplot as plt
from scipy.special import betaincinv  # type: ignore
import hist

myhist = load('hists/stop_2024.scaled')
mcs = myhist['bkg']

signals = ['SMS-2Stop-Par-mStop-1000', 'SMS-2Stop-Par-mStop-1500']

qcd_TvsQCD = mcs['fj1TvsQCD']['QCD Multijet'][{'region': 'cat1_preselection', 'systematic': 'nominal'}].values()
sig_TvsQCD = mcs['fj1TvsQCD'][signals[0]][{'region': 'cat1_preselection', 'systematic': 'nominal'}].values()
second_TvsQCD = mcs['fj1TvsQCD'][signals[1]][{'region': 'cat1_preselection', 'systematic': 'nominal'}].values()
print(len(qcd_TvsQCD))
xaxis = np.linspace(0, 1, 50)
print(len(xaxis))

signal_eff = []
second_eff = []
signal_sum = np.sum(sig_TvsQCD)
second_sum = np.sum(second_TvsQCD)
qcd_rej = []
qcd_sum = np.sum(qcd_TvsQCD)

for cut in xaxis:
    signal_eff.append(np.sum(sig_TvsQCD[xaxis >= cut]) / signal_sum)
    second_eff.append(np.sum(second_TvsQCD[xaxis >= cut]) / second_sum)
    qcd_rej.append(1 - (np.sum(qcd_TvsQCD[xaxis >= cut]) / qcd_sum))

signal_eff.append(0)
second_eff.append(0)
qcd_rej.append(1)

## calculate AUC
auc = np.round(np.trapz(signal_eff, qcd_rej),3)
print("AUC for TvsQCD:", auc)
auc_second = np.round(np.trapz(second_eff, qcd_rej),3)
print("AUC for second signal TvsQCD:", auc_second)

plt.style.use(mplhep.style.CMS)
fig, ax = plt.subplots(figsize=(8, 8))
mplhep.cms.label(llabel='Simulation', rlabel='(13.6 TeV)',fontsize=30, ax=ax)
ax.plot(qcd_rej, signal_eff, label=r'$m_{\tilde{t}}$ = 1000 GeV'+"\n"+r'$m_{\tilde{\chi}}$ = 300 GeV'+"\n"+"AUC = "+str(auc), color='red', linewidth=3.5)
ax.plot(qcd_rej, second_eff, label=r'$m_{\tilde{t}}$ = 1500 GeV'+"\n"+r'$m_{\tilde{\chi}}$ = 100 GeV'+"\n"+"AUC = "+str(auc_second), color='darkred', linewidth=3.5, linestyle='--')
ax.set_ylabel('Signal Efficiency')
ax.set_xlabel('QCD Background Rejection')
ax.set_xlim(0, 1.05)
ax.set_ylim(0, 1.05)
#ax.set_title('Signal Efficiency vs QCD Rejection using TvsQCD')
ax.legend(fontsize=24, loc='lower left')
ax.grid()
fig.savefig('stop_TvsQCD_ROC_vsQCD.png')

# against Zto2Nu

zto2nu_TvsQCD = mcs['fj1TvsQCD']['Z (inv)'][{'region': 'cat1_preselection', 'systematic': 'nominal'}].values()
zto2nu_rej = []
zto2nu_sum = np.sum(zto2nu_TvsQCD)
for cut in xaxis:
    zto2nu_rej.append(1 - (np.sum(zto2nu_TvsQCD[xaxis >= cut]) / zto2nu_sum))
zto2nu_rej.append(1)

# calculate AUC
auc_zto2nu = np.round(np.trapz(signal_eff, zto2nu_rej),3)
print("AUC for TvsQCD vs Zto2Nu:", auc_zto2nu)
auc_second_zto2nu = np.round(np.trapz(second_eff, zto2nu_rej),3)
print("AUC for second signal TvsQCD vs Zto2Nu:", auc_second_zto2nu)

fig, ax = plt.subplots(figsize=(8, 8))
mplhep.cms.label(llabel='Simulation', rlabel='(13.6 TeV)',fontsize=30, ax=ax)
ax.plot(zto2nu_rej, signal_eff, label=r'$m_{\tilde{t}}$ = 1000 GeV'+"\n"+r'$m_{\tilde{\chi}}$ = 300 GeV'+"\n"+"AUC = "+str(auc_zto2nu), color='red', linewidth=3.5)
ax.plot(zto2nu_rej, second_eff, label=r'$m_{\tilde{t}}$ = 1500 GeV'+"\n"+r'$m_{\tilde{\chi}}$ = 100 GeV'+"\n"+"AUC = "+str(auc_second_zto2nu), color='darkred', linewidth=3.5, linestyle='--')
ax.set_ylabel('Signal Efficiency')
ax.set_xlabel('Z→νν Background Rejection')
ax.set_xlim(0, 1.05)
ax.set_ylim(0, 1.05)
#ax.set_title('Signal Efficiency vs Zto2Nu Rejection using TvsQCD')
ax.legend(fontsize=24, loc='lower left')
ax.grid()
fig.savefig('stop_TvsQCD_ROC_vsZto2Nu.png')

# against WtoLNu
wtoLnu_TvsQCD = mcs['fj1TvsQCD']['W (lnu)'][{'region': 'cat1_preselection', 'systematic': 'nominal'}].values()
wtoLnu_rej = []
wtoLnu_sum = np.sum(wtoLnu_TvsQCD)
for cut in xaxis:
    wtoLnu_rej.append(1 - (np.sum(wtoLnu_TvsQCD[xaxis >= cut]) / wtoLnu_sum))
wtoLnu_rej.append(1)
# calculate AUC
auc_wtoLnu = np.round(np.trapz(signal_eff, wtoLnu_rej),3)
print("AUC for TvsQCD vs WtoLNu:", auc_wtoLnu)
auc_second_wtoLnu = np.round(np.trapz(second_eff, wtoLnu_rej),3)
print("AUC for second signal TvsQCD vs WtoLNu:", auc_second_wtoLnu)
fig, ax = plt.subplots(figsize=(8, 8))
mplhep.cms.label(llabel='Simulation', rlabel='(13.6 TeV)',fontsize=30, ax=ax)
ax.plot(wtoLnu_rej, signal_eff, label=r'$m_{\tilde{t}}$ = 1000 GeV'+"\n"+r'$m_{\tilde{\chi}}$ = 300 GeV'+"\n"+"AUC = "+str(auc_wtoLnu), color='red', linewidth=3.5)
ax.plot(wtoLnu_rej, second_eff, label=r'$m_{\tilde{t}}$ = 1500 GeV'+"\n"+r'$m_{\tilde{\chi}}$ = 100 GeV'+"\n"+"AUC = "+str(auc_second_wtoLnu), color='darkred', linewidth=3.5, linestyle='--')
ax.set_ylabel('Signal Efficiency')
ax.set_xlabel('W→lν Background Rejection')
ax.set_xlim(0, 1.05)
ax.set_ylim(0, 1.05)
#ax.set_title('Signal Efficiency vs WtoLNu Rejection using TvsQCD')
ax.legend(fontsize=24, loc='lower left')
ax.grid()
fig.savefig('stop_TvsQCD_ROC_vsWtoLNu.png')