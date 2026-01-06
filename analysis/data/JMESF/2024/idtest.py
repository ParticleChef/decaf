from correctionlib import _core
import awkward as ak

# jet energy correction levels
L1FastJet = "Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi"
L2Relative = "Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi"

## Testing inputs
Jet_area = ak.Array([0.5, 0.7, 0.3])
Jet_eta = ak.Array([0.2, -1.5, 2.3])
Jet_pt = ak.Array([50.0, 100.0, 30.0])
rho = ak.Array([20.0, 15.0, 10.0])

# Load the correction objects
correction_L1FastJet = _core.CorrectionSet.from_file("jet_jerc.json.gz")