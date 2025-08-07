import os
import sys
import argparse

parser = argparse.ArgumentParser(description='MapReduce')
parser.add_argument('-i', '--input', type=str, help='input file')
args = parser.parse_args()

# remove the *.reduced, *.merged and *.scaled files
os.system("rm -f "+ args.input + "/*.reduced")
os.system("rm -f "+ args.input + "/*.merged")


os.system("python3 reduce.py -f "+ args.input)
os.system("python3 merge.py -f "+ args.input)
os.system("python3 scale.py -f "+ args.input + ".merged")
