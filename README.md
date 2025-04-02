<img src="https://user-images.githubusercontent.com/10731328/193421563-cf992d8b-8e5e-4530-9179-7dbd507d2e02.png" width="350"/>

# **D**ark matter **E**xperience with the **C**offea **A**nalysis **F**ramework

---

## Initial Setup

First, log into an LPC node:

```
ssh -L 9094:localhost:9094 <USERNAME>@cmslpc-sl7.fnal.gov
```

The command will also start forwarding the port 9094 (or whatever number you choose)to be able to use applications like jupyter once on the cluster. Then move into your `nobackup` area on `uscms_data`:

```
cd /uscms_data/d?/<USERNAME>
```

where '?' can be [1,2,3]. Install `CMSSW_11_3_4` (Note: `CMSSW 11_3_X` runs on slc7, which can be setup using apptainer on non-slc7 nodes ([see detailed instructions](https://cms-sw.github.io/singularity.html)):

```
#cmssw-el7 # uncomment this line if not on an slc7 node
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
```

Install `combine` ([see detailed instructions](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#installation-instructions)):

```
cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.1.0 # current recommeneded tag (Jan 2024)
scramv1 b clean; scramv1 b # always make a clean build
```

Fork this repo on github and clone it into your `CMSSW_11_3_4/src` directory:

```
cd $CMSSW_BASE/src
git clone https://github.com/<USERNAME>/decaf.git
cd decaf
git switch UL
```

If you want to change run3 branch, following this command
```
git switch run3
```

Then, setup the proper dependences:

```
source setup.sh
```

This script installs the necessary packages as user packages (Note: Pip gave errors when running `setup.sh` for the first time, but it seemed to install everything just fine. No errors showed up when running `setup.sh` a second time.). This is a one-time setup. When you log in next just do:

```
#cmssw-el7 # uncomment this line if not on an slc7 node
cd CMSSW_11_3_4/src
cmsenv
cd decaf
source env.sh
```

By running this script you will also initialize your grid certificate (Note: `setup.sh` also runs `env.sh`). This requires you to save your grid certificate password in `$HOME/private/$USER.txt`. Alternatively, you can comment this out and initialize it manually every time.

---

## Listing Input Files

The list of input files for the analyzer can be generated as a JSON file using the `macros/list.py` script. This script will run over the datasets listed in `data/process.py`, find the list of files for each dataset, “pack” them into small groups for condor jobs, and output the list of groups as a JSON file in `metadata/`.

The options for this script are:

- `-d` (`--dataset`)

Select a specific dataset to pack. By default, it will run over all datasets in `process.py`.

- `-y` (`--year`)

Data year. Options are `2016pre`, `2016post`, `2017`, and `2018`.

- `-m` (`--metadata`)

Name of metadata output file. Output will be saved in `metadata/<NAME>.json`

- `-p` (`--pack`)

Size of file groups. The smaller the number, the more condor jobs will run. The larger the number, the longer each condor job will take. We tend to pick `32`, but the decision is mostly arbitrary.

- `-s` (`--special`)

Size of file groups for special datasets. For a specific dataset, use a different size with respect to the one established with `--pack`. The syntax is `-s <DATASET>:<NUMBER>`.

- `-c` (`--custom`)

Boolean to decide to use public central NanoAODs (if `False`) or private custom NanoAODs (if `True`). Default is `False`.

As an example, to generate the JSON file for all 2017 data:

```
python3 macros/list.py -y 2017 -m 2017 -p 32
```

As a reminder, this script assumes that you are in the `decaf/analysis` directory when running. The output above will be saved in `metadata/2017.json`.

If using the `--custom` option, the script can take several hours to run, so it’s best to use a process manager such as `nohup` or `tmux` to avoid the program crashing in case of a lost connection. For example

```
nohup python3 macros/list.py -y 2017 -m 2017 -p 32 -c &
```

The `&` option at the end of the command lets it run in the background, and the std output and error is saved in `nohup.out`. 

The `nohup` command is useful and recommended for running most scripts, but you may also use tools like `tmux` or `screen`.

---

## Computing MC b-Tagging Efficiencies

MC b-tagging efficiencies are needed by most of the analyses to compute the b-tag event weight, once such efficiencies are corrected with the POG-provided b-tag SFs. To compute them, we first need to run the `common` module in `util`:

```
python utils/common.py
```

This will generate a series of auxiliary functions and information, like the AK4 b-tagging working points, and it will save such information in a `.coffea` file in the `data` folder. AK4 b-tagging working points are essential to measure the MC efficiencies and they are used by the `btag` processor in the `processors` folder. To generate the processor file: 

```
python3 processors/btageff.py -y 2018 -m 2018 -n 2018
```

The options for this script are:

- `-y` (`--year`)

Data year. Options are `2016pre`, `2016post`, `2017`, and `2018`.

- `-m` (`--metadata`)

Metadata file to be used in input.

- `-n` (`--name)

Name of the output processor file. In this case, it will generate a file called `btageff2018.processor` stored in the `data` folder.


To run the processor:

```
python3 run.py -p btageff2018 -m 2018 -d QCD
```

With this command you will run the `btag2018` processor over QCD MC datasets as defined by the `2018` metadata file. You will see a printout like:

Processing: QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8____4_
  Preprocessing 100% ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 32/32 [ 0:01:28 < 0:00:00 | ?   file/s ]
Merging (local) 100% ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 31/31 [ 0:00:23 < 0:00:00 | ? merges/s ]

This means an output file with histograms as defined in the btag processor file has been generated. In this case a folder called `btageff2018` inside the `hists` folder has been created. Inside this folder you can see a file called `QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8____4_.futures`, that stores the histograms. To take advantage of the parallelism offered by the HTCondor job scheduler, the `run_condor.py` script can be used:

```
python3 run_condor.py -p btag2018 -m 2018 -d QCD -c kisti -t -x
```

The options for this script are the same as for `run.py`, with the addition of:

- `-c` (`--cluster`)

Specifies which cluster you are using.  At the moments supports `lpc` or `kisti`.

- `-t` (`--tar`)
  
Tars the local python environment and the local CMSSW folder. 

- `-x` (`--copy`)

Copies these two tarballs to your EOS area. For example, to run the same setup but for a different year you won’t need to tar and copy again. You can simply do: `python run_condor.py -p btag2017 -m 2017 -d QCD -c kisti`

You can check the status of your HTCondor jobs by doing:

```
condor_q <YOUR_USERNAME>
```

After obtaining all the histograms, a first step of data reduction is needed. This step is achieved by running the `reduce.py` script:

```
python reduce.py -f hists/btag2018
```

The options of this script are:

TO BE LISTED

All the different datasets produced at the previous step will be reduced. A different file for each variable for each reduced dataset will be produced. For example, the command above will produce the following reduced files:

```
hists/btageff2018/deepcsv--QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_120to170_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_15to30_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_170to300_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_300to470_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_30to50_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_470to600_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_50to80_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_600to800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_800to1000_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepcsv--QCD_Pt_80to120_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_120to170_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_15to30_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_170to300_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_300to470_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_30to50_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_470to600_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_50to80_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_600to800_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_800to1000_TuneCP5_13TeV_pythia8.reduced
hists/btageff2018/deepflav--QCD_Pt_80to120_TuneCP5_13TeV_pythia8.reduced
```

This step can be run in HTCondor by using the `reduce_condor.py` script. The `reduce_condor.py` script has the same options of `reduce.py`, with addition of the same `--cluster`, `--tar`, and `--copy` options descibed above when discussing `run_condor.py`.

A second step of data reduction is needed to merge all the `.reduced` files corresponding to a single variable. This is achieved by using the `merge.py` script:

```
python3 merge.py -f hists/btageff2018
```

The options of this script are:

TO BE LISTED

This command will produce the following files:

```
hists/btageff2018/deepcsv.merged  hists/btageff2018/deepflav.merged
```

The same script can be used to merge the the files corresponding to each single variable into a single file, using the `-p` or `—postprocess` option:

```
python3 merge.py -f hists/btageff2018 -p
```

Also this step can be run in HTCondor by using the `merge_condor.py` script. The `merge_condor.py` script has the same options of `merge.py`, with addition of the same `--cluster`, `--tar`, and `--copy` options descibed above when discussing `run_condor.py`.

---

This README is a work in progress
