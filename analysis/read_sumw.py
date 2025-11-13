import uproot
import numpy as np
import awkward as ak
from coffea.util import load, save
import mplhep
import matplotlib.pyplot as plt
import hist
from scipy.special import betaincinv # type: ignore

myhist = load('hists/stop_2024.merged')

tot = myhist['sumw']
for key in tot.keys():
    print(key, tot[key])