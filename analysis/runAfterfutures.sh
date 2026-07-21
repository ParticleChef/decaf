#!/bin/bash

pname=${1}

python3 reduce_new.py -f hists/${1}
python3 merge.py -f hists/${1}
python3 macros/scale.py -f hists/${1}.merged

