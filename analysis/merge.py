import cloudpickle
import pickle
import gzip
import os
import numpy as np
from coffea.util import load, save
from helpers.futures_patch import patch_mp_connection_bpo_17560

def merge(folder,variable=None, exclude=None):

     lists = {}
     for filename in os.listdir(folder):
          if '.reduced' not in filename: continue
          if filename.split('--')[0] not in lists: lists[filename.split('--')[0]] = []
          lists[filename.split('--')[0]].append(folder+'/'+filename)

     for var in lists.keys():
          tmp={}
          if variable is not None:
               if not any(v==var for v in variable.split(',')): continue
          if exclude is not None:
               if any(v==var for v in exclude.split(',')): continue
          print(lists[var])
          for filename in lists[var]:
               print('Opening:',filename)
               hin = load(filename)
               if var not in tmp: tmp[var]={}
               if filename.split('--')[1] not in tmp[var]: tmp[var][filename.split('--')[1].replace('.reduced','')]=hin[var]
               del hin
          print(tmp)
          save(tmp, folder+'/'+var+'.merged')


def postprocess(folder):
     
     variables = []
     for filename in os.listdir(folder):
          if '.merged' not in filename: continue
          #if '--' not in filename: continue
          if filename.split('.')[0] not in variables: variables.append(filename.split('.')[0])
     variables = list(filter(None, variables))

     hists = {}
     for variable in variables:
          filename = folder+'/'+variable+'.merged'
          print('Opening:',filename)
          hin = load(filename)
          hists.update(hin)
     print(hists)
     save(hists,folder+'.merged')
     
     

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-v', '--variable', help='variable', dest='variable', default=None)
    parser.add_option('-e', '--exclude', help='exclude', dest='exclude', default=None)
    (options, args) = parser.parse_args()

    patch_mp_connection_bpo_17560()    
    merge(options.folder,options.variable,options.exclude)
    postprocess(options.folder)
