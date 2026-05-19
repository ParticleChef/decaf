#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from coffea.util import load
import hist

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
INPUT = "hists/stop_new2024.scaled"
TARGET_VARIABLE = "metpt"
TARGET_SYSTS = ["pileupUp", "pileupDown"]

REGIONS = [
    "cat1_preselection",
    "cat2_LLCR_highDeltaM",
    "cat3_QCDCR_highDeltaM",
    "cat4_GCR_highDeltaM",
    "cat5_DY2E_highDeltaM",
    "cat6_DY2M_highDeltaM",
    "cat7_SR_highDeltaM",
]

STACKS = [
    "VV",
    "Single Top",
    "TT",
    "DY",
    "Gamma + Jets",
    "W (lnu)",
    "Z (inv)",
    "QCD Multijet",
]

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def reduce_to_1d(vals, vars_=None):
    while vals.ndim > 1:
        vals = vals.sum(axis=0)
        if vars_ is not None:
            vars_ = vars_.sum(axis=0)
    return vals, vars_


def get_rebinned_hist(h, key):
    if key == "j1pt":
        h = h[{"j1pt": hist.rebin(5)}]
    elif key == "j2pt":
        h = h[{"j2pt": hist.rebin(5)}]
    elif key == "fj1pt":
        h = h[{"fj1pt": hist.rebin(5)}]
    elif key == "fj1TvsQCD":
        h = h[{"fj1TvsQCD": hist.rebin(10)}]
    elif key == "metpt_10GeVbins":
        h = h[{"metpt_10GeVbins": hist.rebin(5)}]
    elif key == "recoilpt":
        h = h[{"recoilpt": hist.rebin(5)}]
    return h


def get_hist_safe(container, key, process, region, systematic):
    try:
        h = container[key][process][{"region": region, "systematic": systematic}]
        h = get_rebinned_hist(h, key)
        return h
    except Exception:
        return None


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
myhist = load(INPUT)
bkg = myhist["bkg"]

print(f"[INFO] variable = {TARGET_VARIABLE}")
print(f"[INFO] systematics = {TARGET_SYSTS}")
print(f"[INFO] regions = {REGIONS}")
print()

all_region_results = {}

for region in REGIONS:
    print("\n" + "#" * 120)
    print(f"[REGION] {region}")
    print("#" * 120)

    proc_results = []

    for proc in STACKS:
        if proc not in bkg[TARGET_VARIABLE]:
            print(f"[SKIP] {proc}: not found in bkg[{TARGET_VARIABLE}]")
            continue

        h_nom = get_hist_safe(bkg, TARGET_VARIABLE, proc, region, "nominal")
        if h_nom is None:
            print(f"[SKIP] {proc}: nominal hist missing in {region}")
            continue

        bins = h_nom.axes[-1].edges
        vals_nom, _ = reduce_to_1d(h_nom.values(), None)
        yield_nom = np.sum(vals_nom)

        print("=" * 100)
        print(f"[PROCESS] {proc}")
        print(f"  nominal yield = {yield_nom:.6f}")

        proc_dict = {
            "process": proc,
            "nominal": yield_nom,
        }

        for syst in TARGET_SYSTS:
            h_var = get_hist_safe(bkg, TARGET_VARIABLE, proc, region, syst)
            if h_var is None:
                print(f"  [MISS] {syst}: hist missing")
                proc_dict[syst] = None
                continue

            vals_var, _ = reduce_to_1d(h_var.values(), None)
            yield_var = np.sum(vals_var)

            rel_yield = 100.0 * (yield_var - yield_nom) / yield_nom if yield_nom != 0 else np.nan
            diff = vals_var - vals_nom

            with np.errstate(divide="ignore", invalid="ignore"):
                rel_bin = np.divide(
                    diff,
                    vals_nom,
                    out=np.zeros_like(diff, dtype=float),
                    where=(vals_nom != 0)
                ) * 100.0

            max_abs_bin_idx = int(np.argmax(np.abs(diff)))
            max_rel_bin_idx = int(np.argmax(np.abs(rel_bin)))

            print(f"  {syst:10s} yield = {yield_var:.6f} ({rel_yield:+.3f}%)")
            print(
                f"    max abs bin change : bin {max_abs_bin_idx:2d} "
                f"[{bins[max_abs_bin_idx]:.1f}, {bins[max_abs_bin_idx+1]:.1f}] "
                f"nom={vals_nom[max_abs_bin_idx]:.6f} "
                f"var={vals_var[max_abs_bin_idx]:.6f} "
                f"diff={diff[max_abs_bin_idx]:+.6f}"
            )
            print(
                f"    max rel bin change : bin {max_rel_bin_idx:2d} "
                f"[{bins[max_rel_bin_idx]:.1f}, {bins[max_rel_bin_idx+1]:.1f}] "
                f"nom={vals_nom[max_rel_bin_idx]:.6f} "
                f"var={vals_var[max_rel_bin_idx]:.6f} "
                f"reldiff={rel_bin[max_rel_bin_idx]:+.3f}%"
            )

            top_idx = np.argsort(np.abs(diff))[-10:][::-1]
            print("    top 10 changed bins:")
            for i in top_idx:
                if abs(diff[i]) < 1e-12:
                    continue
                print(
                    f"      bin {i:2d} [{bins[i]:.1f}, {bins[i+1]:.1f}] : "
                    f"nom={vals_nom[i]:.6f}, var={vals_var[i]:.6f}, "
                    f"diff={diff[i]:+.6f}, rel={rel_bin[i]:+.3f}%"
                )

            proc_dict[syst] = {
                "yield": yield_var,
                "rel_percent": rel_yield,
                "max_abs_bin_idx": max_abs_bin_idx,
                "max_abs_diff": float(diff[max_abs_bin_idx]),
                "max_rel_bin_idx": max_rel_bin_idx,
                "max_rel_diff_percent": float(rel_bin[max_rel_bin_idx]),
            }

        proc_results.append(proc_dict)

    all_region_results[region] = proc_results

    print("\n" + "-" * 100)
    print(f"[SUMMARY: {region}] worst offenders by yield change")
    for syst in TARGET_SYSTS:
        valid = []
        for x in proc_results:
            if isinstance(x.get(syst), dict):
                valid.append((x["process"], abs(x[syst]["rel_percent"]), x[syst]["rel_percent"]))
        valid = sorted(valid, key=lambda t: t[1], reverse=True)

        print(f"\n  systematic = {syst}")
        if not valid:
            print("    no valid processes")
            continue

        for proc, absrel, signedrel in valid:
            print(f"    {proc:15s} {signedrel:+10.3f}%")

    print("\n" + "-" * 100)
    print(f"[SUMMARY: {region}] worst offenders by max bin relative change")
    for syst in TARGET_SYSTS:
        valid = []
        for x in proc_results:
            if isinstance(x.get(syst), dict):
                valid.append((x["process"], abs(x[syst]["max_rel_diff_percent"]), x[syst]["max_rel_diff_percent"]))
        valid = sorted(valid, key=lambda t: t[1], reverse=True)

        print(f"\n  systematic = {syst}")
        if not valid:
            print("    no valid processes")
            continue

        for proc, absrel, signedrel in valid:
            print(f"    {proc:15s} {signedrel:+10.3f}%")

print("\n" + "#" * 120)
print("[GLOBAL SUMMARY] worst yield change across all regions")
print("#" * 120)

for syst in TARGET_SYSTS:
    global_valid = []
    for region, proc_results in all_region_results.items():
        for x in proc_results:
            if isinstance(x.get(syst), dict):
                global_valid.append(
                    (region, x["process"], abs(x[syst]["rel_percent"]), x[syst]["rel_percent"])
                )

    global_valid = sorted(global_valid, key=lambda t: t[2], reverse=True)

    print(f"\n  systematic = {syst}")
    for region, proc, absrel, signedrel in global_valid[:20]:
        print(f"    {region:25s} {proc:15s} {signedrel:+10.3f}%")

print("\n" + "#" * 120)
print("[GLOBAL SUMMARY] worst max-bin relative change across all regions")
print("#" * 120)

for syst in TARGET_SYSTS:
    global_valid = []
    for region, proc_results in all_region_results.items():
        for x in proc_results:
            if isinstance(x.get(syst), dict):
                global_valid.append(
                    (
                        region,
                        x["process"],
                        abs(x[syst]["max_rel_diff_percent"]),
                        x[syst]["max_rel_diff_percent"]
                    )
                )

    global_valid = sorted(global_valid, key=lambda t: t[2], reverse=True)

    print(f"\n  systematic = {syst}")
    for region, proc, absrel, signedrel in global_valid[:20]:
        print(f"    {region:25s} {proc:15s} {signedrel:+10.3f}%")