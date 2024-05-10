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

## read merged file
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="input file")
parser.add_argument("-y", type=str, help="year", default="2018", choices=["2016", "2017", "2018"])
parser.add_argument("-b", type=str, help="datablind", default="", choices=["", "blind"])
args = parser.parse_args()
class Mymaker:
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
        self.outdir = self.workingdir + "/datacards/" + self.fname
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
    def makecard(self):
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
        regions = ['sr']
        stacks = [r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV','Single Top', 'Z ($\ell\ell$) + Jets',r'Z ($\nu\nu$) + Jets'] #'$\gamma$ + Jets']
        legs = ['WJets','TT','QCD','VV','ST','ZJets','Zinv','GJets']
        for sign in sig['fj1pt'].keys():
            cardname = self.outdir + f"/datacard_{sign}_{year}.txt" 
            card = open(cardname, "w")
            card.write("imax 1\n")
            card.write(f"jmax {len(stacks)}\n")
            card.write("kmax 1\n") #################### WE NEED TO CHANGE!!!!!!!!
            card.write("---------------\n")
            ### This line will be filled for shape uncertainties
            ### Maybe Next time
            card.write("---------------\n")
            card.write("bin\t\t")
            for region in regions:
                card.write("{0}\t\t".format(region))
            if blind == "blind":
                card.write("\nobservation -1\n")
            else:
                print("You are trying to unblind the data. Exiting...")
                sys.exit(1)
            card.write("---------------\n")
            card.write("bin\t\t")
            for region in regions:
                for stack in stacks:
                    card.write("{0}\t\t".format(region))
            card.write("sr\n")
            card.write("process\t\t")
            for region in regions:
                for i, stack in enumerate(stacks):
                    card.write("{0}\t\t".format(legs[i]))
            card.write("signal\n")
            card.write("process\t\t")
            for region in regions:
                for i in range(len(stacks)):
                    card.write("{0}\t\t".format(i+1))
            card.write("0\n")
            card.write("rate\t\t")
            for region in regions:
                for i, stack in enumerate(stacks):
                    targ = legs[i],np.sum(bkg['fj1pt'][stack][{'region': region, 'TvsQCD': sum}].values())
                    card.write("{0}\t\t".format(np.sum(bkg['fj1pt'][stack][{'region': region, 'TvsQCD': slice(0.26j, 1.1j, sum)}].values())))
            card.write("{0}\n".format(np.sum(sig['fj1pt'][sign][{'region': 'sr', 'TvsQCD': slice(0.26j, 1.1j, sum)}].values())))
            #card.write("10000\n")
            card.write("----------------------\n")
            card.write("lumi\t\tlnN\t\t")
            for region in regions:
                for stack in stacks:
                    card.write("1.025\t\t")
            card.write("1.025\n")

            card.close()

            
#
if __name__ == "__main__":
    maker = Mymaker(args.i, args.y, args.b)
    maker.makecard()
    print("Datacard is created.")
    print("Exiting...")
    sys.exit(0)
