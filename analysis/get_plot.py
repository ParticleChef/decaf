import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import mplhep
import os, sys

# get limit.csv
df = pd.read_csv("limit.csv")
df = df.sort_values(by="Mediator")

# plot the limits
plt.figure(figsize=(9,9))
plt.style.use(hep.style.CMS)
hep.cms.label(llabel="Work in Progress", rlabel="59.8 fb$^{-1}$ (13 TeV)")
# All values are log10

plt.plot(df["Mediator"], np.log10(df["exp"]), label="Median Expected", color="black", linewidth=2, linestyle="--")
plt.fill_between(df["Mediator"], np.log10(df["m2"]), np.log10(df["p2"]), color="#F5BB54", alpha=1, label="95% expected")
plt.fill_between(df["Mediator"], np.log10(df["m1"]), np.log10(df["p1"]), color="#607641", alpha=1, label="68% expected")

plt.axhline(np.log10(1), color="red", linestyle=":", label="$\mu$ = 1 Exclusion")
plt.text(x=0.05, y=0.9, s = "Dirac, AxialVector", ha='left', va='center', transform=plt.gca().transAxes,fontsize=25, fontweight='bold')
plt.text(x=0.07, y=0.83, s = "$g_{q} = 0.25$, $g_{\chi} = 1.0$", ha='left', va='center', transform=plt.gca().transAxes, fontsize=18)
plt.legend(loc="lower right")


plt.grid()
plt.xlim(500,3000)
plt.xlabel("m$_{med}$ (GeV)")
plt.ylabel("log$_{10}$($\mu$)")

plt.savefig("limit2018.png")
