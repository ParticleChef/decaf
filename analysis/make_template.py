#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import uproot
import hist
from coffea.util import load

# ============================================================
# Config
# ============================================================
INPUT = "hists/stop_new2024.scaled"
OUTPUT = "templates_metpt.root"

VARIABLE = "metpt"

# 저장할 region들
REGIONS = [
    "cat2_LLCR_highDeltaM",
    "cat3_QCDCR_highDeltaM",
    "cat4_GCR_highDeltaM",
    "cat5_DY2E_highDeltaM",
    "cat6_DY2M_highDeltaM",
    "cat7_SR_highDeltaM",
]

# background / signal process 이름
BKG_PROCESSES = [
    "VV",
    "Single Top",
    "TT",
    "DY",
    "Gamma + Jets",
    "W (lnu)",
    "Z (inv)",
    "QCD Multijet",
]

SIG_PROCESSES = [
    "SMS-2Stop-Par-mStop-1000",
    "SMS-2Stop-Par-mStop-1500",
    "SMS-2Stop-Par-mStop-600",
]

# shape systematic 목록
SHAPE_SYSTS = {
    "pileup": ("pileupUp", "pileupDown"),
}

# Combine-friendly process names
PROCESS_NAME_MAP = {
    "VV": "VV",
    "Single Top": "SingleTop",
    "TT": "TT",
    "DY": "DY",
    "Gamma + Jets": "GammaJets",
    "W (lnu)": "Wlnu",
    "Z (inv)": "Zinv",
    "QCD Multijet": "QCD",
    "SMS-2Stop-Par-mStop-1000": "SMS_2Stop_mStop1000",
    "SMS-2Stop-Par-mStop-1500": "SMS_2Stop_mStop1500",
    "SMS-2Stop-Par-mStop-600": "SMS_2Stop_mStop600",
}

# ============================================================
# Helpers
# ============================================================
def get_data_process(region):
    if "GCR" in region or "DY2E" in region:
        return "EGamma"
    elif "DY2M" in region:
        return "Muon"
    else:
        return "JetMET"


def sanitize_name(name):
    if name in PROCESS_NAME_MAP:
        return PROCESS_NAME_MAP[name]
    return (
        name.replace(" ", "")
            .replace("(", "")
            .replace(")", "")
            .replace("+", "")
            .replace("-", "_")
            .replace("/", "_")
    )


def get_hist_safe(container, variable, process, region, systematic):
    try:
        h = container[variable][process][{"region": region, "systematic": systematic}]
        return h
    except Exception:
        return None


def ensure_1d_hist(h):
    """
    region/systematic/process를 slice한 뒤에는 보통 1D hist여야 함.
    혹시 추가 축이 남아 있으면 마지막 축만 남기고 나머지는 sum.
    """
    if h is None:
        return None

    while len(h.axes) > 1:
        # 첫 축을 모두 합치기
        h = h[{h.axes[0].name: sum}]

    return h


def hist_has_content(h):
    if h is None:
        return False
    vals = h.values()
    if vals is None:
        return False
    return np.any(np.isfinite(vals))


def make_empty_like(nominal_hist):
    """
    nominal hist와 같은 binning의 빈 1D weighted hist 생성
    """
    ax = nominal_hist.axes[0]
    if isinstance(ax, hist.axis.Variable):
        new_h = hist.Hist(
            hist.axis.Variable(ax.edges, name=ax.name, label=ax.label),
            storage=hist.storage.Weight()
        )
    elif isinstance(ax, hist.axis.Regular):
        new_h = hist.Hist(
            hist.axis.Regular(ax.size, ax.edges[0], ax.edges[-1], name=ax.name, label=ax.label),
            storage=hist.storage.Weight()
        )
    else:
        raise RuntimeError(f"Unsupported axis type for template writing: {type(ax)}")
    return new_h


def set_negative_bins_to_zero(h):
    """
    Combine용으로 음수 bin content가 있으면 0으로 자름.
    variance는 그대로 두거나, content가 0이면 variance도 0으로 둠.
    """
    vals = np.array(h.values(), dtype=float)
    vars_ = h.variances()
    if vars_ is None:
        vars_ = np.zeros_like(vals, dtype=float)
    else:
        vars_ = np.array(vars_, dtype=float)

    mask = vals < 0
    vals[mask] = 0.0
    vars_[mask] = 0.0

    out = make_empty_like(h)
    view = out.view(flow=False)
    view.value = vals
    view.variance = vars_
    return out


def prepare_template_hist(h, clip_negative=True):
    h = ensure_1d_hist(h)
    if h is None:
        return None
    if clip_negative:
        h = set_negative_bins_to_zero(h)
    return h


# ============================================================
# Main
# ============================================================
def main():
    if not os.path.exists(INPUT):
        raise FileNotFoundError(f"Input file not found: {INPUT}")

    myhist = load(INPUT)
    bkg = myhist["bkg"]
    data = myhist["data"]
    sig = myhist["sig"] if "sig" in myhist else myhist["bkg"]

    all_mc_processes = BKG_PROCESSES + SIG_PROCESSES

    with uproot.recreate(OUTPUT) as fout:
        for region in REGIONS:
            print(f"[INFO] Writing region: {region}")

            # --------------------------------------------------------
            # data_obs
            # --------------------------------------------------------
            data_proc = get_data_process(region)
            h_data = get_hist_safe(data, VARIABLE, data_proc, region, "nominal")
            h_data = prepare_template_hist(h_data, clip_negative=False)

            if h_data is not None:
                fout[f"{region}/data_obs"] = h_data
                print(f"  wrote {region}/data_obs")
            else:
                print(f"  [WARN] missing data for region={region}, process={data_proc}")

            # --------------------------------------------------------
            # MC nominal
            # --------------------------------------------------------
            nominal_cache = {}

            for proc in all_mc_processes:
                source = bkg if proc in bkg.get(VARIABLE, {}) else sig
                h_nom = get_hist_safe(source, VARIABLE, proc, region, "nominal")
                h_nom = prepare_template_hist(h_nom, clip_negative=True)

                if h_nom is None:
                    print(f"  [WARN] missing nominal: region={region}, process={proc}")
                    continue

                proc_name = sanitize_name(proc)
                nominal_cache[proc] = h_nom
                fout[f"{region}/{proc_name}"] = h_nom
                print(f"  wrote {region}/{proc_name}")

            # --------------------------------------------------------
            # Shape systematics
            # --------------------------------------------------------
            for proc in all_mc_processes:
                if proc not in nominal_cache:
                    continue

                source = bkg if proc in bkg.get(VARIABLE, {}) else sig
                proc_name = sanitize_name(proc)
                h_nom = nominal_cache[proc]

                for syst_name, (syst_up, syst_down) in SHAPE_SYSTS.items():
                    h_up = get_hist_safe(source, VARIABLE, proc, region, syst_up)
                    h_down = get_hist_safe(source, VARIABLE, proc, region, syst_down)

                    h_up = prepare_template_hist(h_up, clip_negative=True)
                    h_down = prepare_template_hist(h_down, clip_negative=True)

                    # 없으면 nominal 복사
                    if h_up is None:
                        h_up = h_nom
                    if h_down is None:
                        h_down = h_nom

                    fout[f"{region}/{proc_name}_{syst_name}Up"] = h_up
                    fout[f"{region}/{proc_name}_{syst_name}Down"] = h_down
                    print(f"  wrote {region}/{proc_name}_{syst_name}Up")
                    print(f"  wrote {region}/{proc_name}_{syst_name}Down")

    print(f"\n[INFO] Done. Output written to: {OUTPUT}")
    print("\nExample datacard shapes line:")
    print(f"shapes * * {OUTPUT} $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC")


if __name__ == "__main__":
    main()