import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import correctionlib._core as core
import os

# Load your ROOT file as a NanoEvents object
events = NanoEventsFactory.from_root(
    "/data/mc/Run3Summer22NanoAODv12/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/f7267ea1-1288-4112-a06a-849bf1f36dfa.root",  # <-- change to your file
    schemaclass=NanoAODSchema,
    entry_start=0,
    entry_stop=None
).events()

# Access jets (change 'Jet' to 'FatJet' for AK8, etc.)
jets = events.Jet

# Load the JEC JSON file
jec_json_path = os.path.join(os.path.dirname(__file__), "data/JetMETCorr/2022pre/jet_jerc.json.gz")
cset = core.CorrectionSet.from_file(jec_json_path)

# Set your JEC key (adjust as needed)
jec = "Summer19UL16_V7_MC"
lvl = "L2Relative"  # or "L1L2L3Res", etc.
algo = "AK4PFchs"
key = f"{jec}_{lvl}_{algo}"

# Get the correction object
jec_correction = cset[key]

# Prepare the input variables for the correction
# These must match the order in the JSON
inputs = [
    ak.to_numpy(jets.pt),    # JetPt
    ak.to_numpy(jets.eta),   # JetEta
    ak.to_numpy(jets.area),  # JetA
    ak.to_numpy(jets.phi),   # JetPhi
    ak.to_numpy(events.fixedGridRhoFastjetAll)  # Rho
]

# Evaluate the correction factor for each jet
jec_factors = jec_correction.evaluate(*inputs)

# Apply the correction to jet pt and mass
jets["pt_corr"] = jets.pt * jec_factors
jets["mass_corr"] = jets.mass * jec_factors

# Now jets.pt_corr and jets.mass_corr are the corrected values
print(jets.pt_corr)
print(jets.mass_corr)
