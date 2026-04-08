## Changes: 2026-03-10 (First 2024 b-tagging SF combination using kinfit and ptRel methods)

Merge Request: [!3](https://gitlab.cern.ch/cms-analysis-corrections/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/-/merge_requests/3)

The 2024 SF corrections are updated to include the first combination of b-tagging SFs in the 2024 era computed with the kinematic fit method in the ttbar phase space and the ptRel method in the QCD phase space. New corrections with names `UParTAK4_comb`, `UParTAK4_mujets` and `UParTAK4_light` are added to `btagging.json.gz` to include the combined SF.

The description of the corrections is the following:
- `UParTAK4_comb`: combined b-tagging SF for b and c jets (hadronFlavour=4,5) obtained with the combination of the kinfit and ptRel methods
- `UParTAK4_mujets`: b-tagging SF for b and c jets (hadronFlavour=4,5) obtained with the ptRel method only
- `UParTAK4_light`: b-tagging SF for light jets (hadronFlavour=0) obtained with the negTag (DY phase space) method only

The preliminary SF file `btagging_preliminary.json.gz` is deleted, since now the combined SF are available.
The format of the SFs provided in `btagging.json.gz` now is the same as the SFs provided in the 22+23 eras.

In the future, the `UParTAK4_comb` SFs will be updated by combining the remaining calibration methods.
Similarly, the `UParTAK4_mujets` SFs will be updated by combining the remaining calibration methods in the QCD phase space.
Finally, also the `UParTAK4_light` SFs will be updated when the SFs obtained with the negTag method in the QCD phase space will be available.
