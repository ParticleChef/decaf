import correctionlib    
def check_json(path):
    year = '2022pre'
    evaluator = correctionlib.CorrectionSet.from_file(path)
    for corr in evaluator.values():
        #if not corr.name == "NUM_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose": continue
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        #print(f"{corr.name}")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")


check_json("/home/twkim/NPS_Stop_Analysis/decaf/analysis/data/JMESF/2022pre/fatJet_jerc.json.gz")     