import uproot
import numpy as np
import awkward as ak
from coffea.util import load, save
import mplhep
import matplotlib.pyplot as plt
import hist

myhist = load("hists/stop_2022pre/EGammaC____1_.futures")

print(myhist)