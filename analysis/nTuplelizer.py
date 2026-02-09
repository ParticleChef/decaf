#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make per-file NanoAOD nTuples keeping ALL branches, selecting ONLY by (OR of) signal triggers.

- Input : /data/data/2024*/JetMET*/*.root
- Output: $(pwd)/ntuple_out/<Era>/ <basename>_trig.root
- Keeps: Events + Runs + LuminosityBlocks (same branches as input)
- ROOT-version friendly: no Experimental progress API, no Report() indexing
- Resume: 이미 생성된 출력(Events entries > 0)은 스킵
"""

import ROOT, glob, os, time, re, fnmatch

ROOT.ROOT.EnableImplicitMT()  # use all cores

# ---------- Config ----------
IN_GLOB  = "/data/data/2024*/JetMET*/*.root"
TREENAME = "Events"

# Output root under current working dir, separated by Era (e.g. 2024A, 2024B, ...)
OUT_ROOT = os.path.join(os.getcwd(), "ntuple_out")
os.makedirs(OUT_ROOT, exist_ok=True)

# OR-list of triggers to pass (edit if needed)
SIGNAL_TRIGGERS = [
    "HLT_PFMET120_PFMHT120_IDTight",
    "HLT_PFMET130_PFMHT130_IDTight",
    "HLT_PFMET140_PFMHT140_IDTight",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
]

# Snapshot options
opts_recreate = ROOT.RDF.RSnapshotOptions()
opts_recreate.fMode = "RECREATE"
opts_recreate.fCompressionAlgorithm = ROOT.kZSTD
opts_recreate.fCompressionLevel = 6
opts_recreate.fAutoFlush = 30_000_000

opts_update = ROOT.RDF.RSnapshotOptions()
opts_update.fMode = "UPDATE"
opts_update.fCompressionAlgorithm = ROOT.kZSTD
opts_update.fCompressionLevel = 6
opts_update.fAutoFlush = 30_000_000
# ----------------------------


def extract_era(path: str) -> str:
    """Extract '2024X' from the file path; return 'Unknown' if not found."""
    m = re.search(r"/(2024[A-Z])/", path)
    return m.group(1) if m else "Unknown"


def build_trigger_expr_for_file(infile: str):
    """Return 'HLT_A || HLT_B || ...' using only triggers that EXIST in this file.
       If none exist, return None (pass-through)."""
    rdf_probe = ROOT.RDataFrame(TREENAME, infile)
    cols = set(rdf_probe.GetColumnNames())
    use = [t for t in SIGNAL_TRIGGERS if t in cols]
    return " || ".join(use) if use else None


def _count_entries_fast(filename: str, treename: str) -> int:
    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        if f: f.Close()
        return 0
    t = f.Get(treename)
    n = int(t.GetEntriesFast()) if t else 0
    f.Close()
    return n


def is_valid_output(outfile: str, treename: str = "Events") -> bool:
    """Resume 체크: 파일 존재 + Events 트리 존재 + entries > 0."""
    if not os.path.isfile(outfile):
        return False
    try:
        return _count_entries_fast(outfile, treename) > 0
    except Exception:
        return False


def snapshot_all_trees_like_input(infile, outfile, trigger_expr_or_none):
    """Write Events (RECREATE) + Runs/LuminosityBlocks (UPDATE) to 'outfile'."""
    # Base dataframe and counts
    rdf0 = ROOT.RDataFrame(TREENAME, infile)
    n_total = int(rdf0.Count().GetValue())

    # Apply trigger OR (if any trigger exists in this file)
    rdf = rdf0.Filter(trigger_expr_or_none, "signal_trigger_any") if trigger_expr_or_none else rdf0
    n_pass = int(rdf.Count().GetValue())

    # 1) write Events first (RECREATE so file is created fresh)
    t0 = time.time()
    rdf.Snapshot(TREENAME, outfile, "", opts_recreate)

    # 2) append auxiliary trees (UPDATE so we don't overwrite Events)
    for aux_tree in ("Runs", "LuminosityBlocks"):
        try:
            rdf_aux = ROOT.RDataFrame(aux_tree, infile)
            # Snapshot will safely write zero-entry trees too, but cheap to check:
            if rdf_aux.Count().GetValue() >= 0:
                rdf_aux.Snapshot(aux_tree, outfile, "", opts_update)
        except Exception:
            # Tree may not exist in some files — ignore
            pass

    dt = time.time() - t0
    kept_rate = (n_pass / n_total * 100.0) if n_total else 0.0
    speed = (n_total / dt) if dt > 0 else 0.0
    print(f"      ⏱ {dt:.1f}s | total {n_total:,} → kept {n_pass:,} ({kept_rate:.2f}%) | {speed:,.0f} evt/s")


def main():
    # Collect inputs, EXCLUDING any already-produced *_trig.root
    files = [f for f in sorted(glob.glob(IN_GLOB)) if not fnmatch.fnmatch(f, "*_trig.root")]
    if not files:
        print(f"[WARN] No files matched: {IN_GLOB}")
        return

    # prescan: 이미 유효한 출력이 있으면 큐에서 제외
    to_run = []
    skipped = 0
    for f in files:
        b = os.path.basename(f)
        era = extract_era(f)
        out_dir = os.path.join(OUT_ROOT, era)
        os.makedirs(out_dir, exist_ok=True)
        name, _ = os.path.splitext(b)
        outfile = os.path.join(out_dir, f"{name}_trig.root")
        if is_valid_output(outfile, TREENAME):
            skipped += 1
        else:
            to_run.append((f, outfile))

    print(f"🔹 Found {len(files)} files (excluding *_trig.root) — prescan skip {skipped}, to process {len(to_run)}")

    if not to_run:
        print("✅ All outputs already exist and look valid. Nothing to do.")
        print(f"📂 Output root: {OUT_ROOT}")
        return

    # process sequentially (RDF 내부는 멀티스레드)
    for i, (f, outfile) in enumerate(to_run, 1):
        b = os.path.basename(f)
        expr = build_trigger_expr_for_file(f)
        if expr:
            print(f"[{i}/{len(to_run)}] {b}: triggers → {expr}")
        else:
            print(f"[{i}/{len(to_run)}] {b}: no trigger branches → pass-through (None→True)")

        # 동시 실행 등으로 출력이 생겼을 수 있으니 한 번 더 체크
        if is_valid_output(outfile, TREENAME):
            print(f"      ↪︎ SKIP_EXISTING: {outfile}")
            continue

        snapshot_all_trees_like_input(f, outfile, expr)

    print("✅ All requested files processed.")
    print(f"📂 Output root: {OUT_ROOT}")


if __name__ == "__main__":
    main()
