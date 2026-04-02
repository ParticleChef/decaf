import uproot
import numpy as np
import awkward as ak
from coffea.nanoevents.methods import vector
from coffea.lumi_tools import LumiMask
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import gzip
import json


lumimask = LumiMask('/home/jhong/run3Monotop/decaf/analysis/data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json')

# Muon ID
def isTightMuon(pt, eta, iso, id):
    mask = (pt > 30) & (abs(eta) < 2.4) & (iso < 0.1) & (id == 1)
    return mask

txtfile = 'Run2022postEE.txt'
recoil_bins = np.linspace(0, 1000, 51)

def process_one_file(args):
    i, nfiles, file = args
    print('Processing file %d/%d: %s' % (i+1, nfiles, file))

    max_retries = 3
    timeout_s = 180  # Timeout 3min
    last_err = None

    for attempt in range(1, max_retries + 1):
        try:
            with uproot.open(file, timeout=timeout_s) as f:
                tree = f['Events']

                muon = tree.arrays(
                    ['Muon_pt', 'Muon_eta', 'Muon_phi', 'Muon_pfRelIso04_all', 'Muon_tightId'],
                    how='zip'
                )['Muon']

                met = tree.arrays(['MET_pt', 'MET_phi'], how='zip')

                reference_trigger = tree['HLT_IsoMu24'].array(library='ak')
                met_trigger = (
                    tree['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'].array(library='ak')
                    | tree['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60'].array(library='ak')
                )

                filters = (
                    tree['Flag_goodVertices'].array(library='ak')
                    & tree['Flag_globalSuperTightHalo2016Filter'].array(library='ak')
                    & tree['Flag_EcalDeadCellTriggerPrimitiveFilter'].array(library='ak')
                    & tree['Flag_BadPFMuonFilter'].array(library='ak')
                    & tree['Flag_BadPFMuonDzFilter'].array(library='ak')
                    & tree['Flag_hfNoisyHitsFilter'].array(library='ak')
                    & tree['Flag_eeBadScFilter'].array(library='ak')
                )

                lumi_mask = lumimask(tree['run'].array(), tree['luminosityBlock'].array())

                muon_tight = muon[isTightMuon(muon.pt, muon.eta, muon.pfRelIso04_all, muon.tightId)]
                muon_T = ak.zip({'r': muon_tight.pt, 'phi': muon_tight.phi},
                                with_name='PolarTwoVector', behavior=vector.behavior)
                met_T = ak.zip({'r': met.MET_pt, 'phi': met.MET_phi},
                               with_name='PolarTwoVector', behavior=vector.behavior)

                nMuon = ak.num(muon)
                met100 = met.MET_pt > 100
                leading_mu = ak.firsts(muon_T)
                recoil = (met_T + leading_mu).r
                # None to zero
                recoil = ak.fill_none(recoil, 0)
                recoil70 = (recoil > 70)

                selection = (met100) & (nMuon == 1) & (reference_trigger) & (filters) & (lumi_mask) & (recoil70)
                recoil = recoil[selection]
                met = met[selection]
                met_trigger = met_trigger[selection]


                recoil_triggered = recoil[met_trigger]

                trig = plt.hist(recoil_triggered, bins=recoil_bins, histtype='step',
                                label='Triggered', color='blue')
                tot = plt.hist(recoil, bins=recoil_bins, histtype='step',
                               label='All', color='red')

                plt.clf()
                return trig[0], tot[0]

        except Exception as e:
            last_err = e
            print(f"[WARN] failed (attempt {attempt}/{max_retries}) for {file}: {e}")
            time.sleep(2 * attempt)

    print(f"[SKIP] giving up on {file} after {max_retries} retries. Last error: {last_err}")
    zeros = np.zeros(len(recoil_bins) - 1, dtype=np.float64)
    return zeros, zeros

if __name__ == "__main__":
    total = []
    triggered = []

    with gzip.open('/home/jhong/run3Monotop/decaf/analysis/metadata/for_recoil_trigger.json.gz', 'rt') as f:
        data = json.load(f)

    files = []
    for key in data:
        files.extend(data[key]["files"])

    files_to_run = files #[:10]
    nfiles = len(files_to_run)
    tasks = [(i, nfiles, file) for i, file in enumerate(files_to_run)]

    max_workers = 8

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(process_one_file, t) for t in tasks]
        for fut in as_completed(futures):
            trig0, tot0 = fut.result()
            triggered.append(trig0)
            total.append(tot0)

    triggered = np.sum(triggered, axis=0)
    total = np.sum(total, axis=0)

    efficiency = np.divide(triggered, total, out=np.zeros_like(triggered, dtype=float), where=(total != 0))
    print('Efficiency:', efficiency)

    hep.style.use('CMS')
    plt.figure(figsize=(8, 8))
    plt.step(recoil_bins[:-1], efficiency, where='post', label='Trigger Efficiency',
             color='green', linewidth=2.5)
    plt.xlabel(r"$U\!\!\!/_{T}\ \mathrm{(GeV)}$")
    plt.ylabel('Trigger Efficiency')
    plt.legend()
    plt.grid()
    plt.savefig('trigger_efficiency_22post_mc.png')

    with uproot.recreate("recoil_trigger_eff_22post_mc.root") as f:
        f["recoil_efficiency"] = (efficiency, recoil_bins)
        f["triggered"] = (triggered, recoil_bins)
        f["total"] = (total, recoil_bins)
