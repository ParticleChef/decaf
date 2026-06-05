<img src="https://user-images.githubusercontent.com/10731328/193421563-cf992d8b-8e5e-4530-9179-7dbd507d2e02.png" width="350"/>

# **D**ark matter **E**xperience with the **C**offea **A**nalysis **F**ramework

---

*README updating now*

## Initial Setup for `coffea 0.7.22` with Python 3.8.x

This branch of decaf relies on an older Python/coffea software stack that is no longer available in the default CMS AlmaLinux environment. To ensure compatibility, users should create a dedicated Python 3.8.18 virtual environment and install the required dependencies as described below.

To avoid `xrootd` error, **conda** also can be used for setup.

Move to your working area:

```bash
cd /path/to/workdir
```

### Create a Python 3.8.18 environment

Download and build Python 3.8.18:

```bash
wget https://www.python.org/ftp/python/3.8.18/Python-3.8.18.tgz
tar -xvf Python-3.8.18.tgz
cd Python-3.8.18

./configure  --prefix=/path/to/python/python-3.8.18  --enable-optimizations

make -j$(nproc)
make install
```

In `--prefix` option, put the path of whwere you want to install python at.


Create a virtual environment:

```bash
cd ..
mkdir envs/
/path/to/python/python-3.8.18/bin/python3.8 -m venv ./envs/p38
```
You can change the directory for virtual environment directory (`envs/`)


Activate it:

```bash
source ./envs/p38/bin/activate
```

Upgrade pip tools:

```bash
pip install --upgrade pip setuptools wheel
```

`setuptools` must have version less than `65.7.0`. 

### Install analysis dependencies

```bash
pip install \
    coffea==0.7.22 \
    awkward==1.10.5 \
    numpy==1.23.5 \
    uproot==4.3.7 \
    vector==1.3.1 \
    hist==2.7.2

pip install https://github.com/mcremone/rhalphalib/archive/master.zip

pip install xxhash
pip install 'correctionlib[convert]'
pip install tabulate
```

### Clone decaf

```bash
git clone -b run3 https://github.com/ParticleChef/decaf.git
cd decaf
```

### Start the environment

For every new session:

```bash
source /path/to/envs/p38/bin/activate
cd decaf/analysis
```
Ignore the `env.sh` file now.

---

## Listing Input Files

*Will be updated this part soon*

The list of input files for the analyzer can be generated as a JSON file using the `macros/list.py` script. This script will run over the datasets listed in `data/process.py`, find the list of files for each dataset, “pack” them into small groups for condor jobs, and output the list of groups as a JSON file in `metadata/`.

The options for this script are:

- `-d` (`--dataset`)

Select a specific dataset to pack. By default, it will run over all datasets in `process.py`.

- `-y` (`--year`)

Data year. Now run2 is ready but if you want to make Run3 files, please add the campaign array.

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

## Start Analysis
Run object definition file:
```
python3 utils/ids.py
```

Run correction evaluation function file. Corretion files are in `data/` directory.:
```
python3 utils/corrections.py
```

Run b-tag working points file:
```
python3 utils/common.py
```

The main Analyzer file is `processors/hadmonotop_run3.py` file. The btag efficiency processor is `btageff_run3.py`. These should be compile and the processor file must be generated before run. The btageff for MCs listed in `decaf/analysis/run3_datasets_XS.csv` file already saved in `hists` directory
To generate the processor file: 

```
python3 processors/hadmonotop_run3.py -y 2022pre -m 2022_private_v1 -n 2022_0605
```

The options for this script are:

- `-y` (`--year`)

Data year. 

- `-m` (`--metadata`)

Metadata file to be used in input.

- `-n` (`--name)

Name of the output processor file. In this case, it will generate a file called `hadmonotop2022_0605.processor` stored in the `data` folder.


### Run the processor

First, for test run with one job in local:

```
python3 run.py -p hadmonotop2022_0605 -m 2022_private_v1 -d QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8____1_
```

To run with `nohup` in local:
```
python3 nohup_job_new.py -p hadmonotop2022_0605 -m 2022_private_v1
```

To run with `condor`, use the coffea image file.
In `run_condor.py`, add:

```python
+SingularityImage = "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base-almalinux8:0.7.22-py3.8"
```

to the condor submit description before submitting jobs.

Run:
```
python3 run_condor.py  -p hadmonotop2022_0605 -m 2022_private_v1 -t -x
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

