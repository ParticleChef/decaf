import correctionlib
from correctionlib import convert
import awkward as ak

# jet energy correction levels
L1FastJet = "Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi"
L2Relative = "Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi"
L3Absolute = "Summer24Prompt24_V1_MC_L3Absolute_AK4PFPuppi"
L2L3Residual = "Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi"
L1L2L3Residual = "Summer24Prompt24_V1_MC_L1L2L3Res_AK4PFPuppi"

## Testing inputs
Jet_area = ak.Array([0.5, 0.7, 0.3])
Jet_eta = ak.Array([0.2, -1.5, 2.3])
Jet_phi = ak.Array([0.1, -2.0, 0.5])
Jet_pt = ak.Array([50.0, 100.0, 30.0])
rho = ak.Array([20.0, 15.0, 10.0])

# Load the correction objects
evaluator = correctionlib.CorrectionSet.from_file("jet_jerc.json.gz")
corr_L1 = evaluator[L1FastJet]
corr_L2 = evaluator[L2Relative]
corr_L3 = evaluator[L3Absolute]
corr_L2L3Residual = evaluator[L2L3Residual]
test = evaluator["123124fdgf"]
corr_L1L2L3Residual = evaluator["Summer24Prompt24_V1_MC_L1L2L3Res_AK4PFPuppi"]
print(evaluator.keys())

# Apply the L1FastJet correction
corrected_pt_L1 = corr_L1.evaluate(Jet_area, Jet_eta, Jet_pt, rho)
corrected_pt_L2 = corr_L2.evaluate(Jet_eta, Jet_phi, corrected_pt_L1)
corrected_pt_L3 = corr_L3.evaluate(Jet_eta, corrected_pt_L2)
corrected_pt_L2L3 = corr_L2L3Residual.evaluate(Jet_eta, corrected_pt_L3)
#corrected_pt_L1L2L3 = corr_L1L2L3Residual.evaluate(Jet_area, Jet_eta, Jet_pt, rho, Jet_phi)
print("Corrected pt after L1FastJet:", corrected_pt_L1)
print("Corrected pt after L2Relative:", corrected_pt_L2)
print("Corrected pt after L3Absolute:", corrected_pt_L3)
print("Corrected pt after L2L3Residual:", corrected_pt_L2L3)
#print("Corrected pt after L1L2L3Residual:", corrected_pt_L1L2L3)


