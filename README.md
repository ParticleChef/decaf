<img src="https://user-images.githubusercontent.com/10731328/193421563-cf992d8b-8e5e-4530-9179-7dbd507d2e02.png" width="350"/>

# **D**ark matter **E**xperience with the **C**offea **A**nalysis **F**ramework
Following instructions are to run the generation of the histograms directly from NanoAOD. 

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
git clone git@github.com:<USERNAME>/decaf.git
cd decaf
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

This README is a work in progress
