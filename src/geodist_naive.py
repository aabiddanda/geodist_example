#!/usr/local/bin/python3

'''
  Naive categorization of geodist categories
'''

from __future__ import print_function
import sys

import click
import numpy as np 
import gzip as gz
import json

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def generate_bins(endpts):
  assert(np.all(np.array(endpts) < 1.0))
  b = 0.0
  y = len(endpts)
  bins = []
  for x in endpts:
      bins.append((b,x))
      b =x
  bins.append((b,1.0))
  return(bins)

def simple_geodist_binning(X, bins=[(0,0),(0,0.01),(0.01,0.05),(0.05,1.0)]):
  """
      Calculate the **simple** version of Geodist labelings
  """
  output = np.zeros(shape = X.shape, dtype=int)
  i = 1
  for b in bins[1:]:
      idx = np.where((X > b[0]) & (X <= b[1]))
      output[idx] = i
      i += 1
  return(output)

"""
  Generate numeric codes 
"""
def generate_codes(X, n_cats=4, n_pops=5):
  bases = n_cats**np.arange(n_pops)
  codes = np.array((X @ bases) + 1, dtype=int)
  return(codes)

"""
  Simply concatenate the categories 
"""
def concat_categories(X):
  Y = X.astype('str').tolist()
  codes = [''.join(row) for row in Y]
  return(codes)

def gen_geodist(infile, outfile='test.geodist.txt.gz', bins="[0.0,0.05]"):
  """ Computing the underlying geographic distribution """
  bin_endpts = list(map(float, bins.strip('[]').split(',')))
  bins = generate_bins(bin_endpts)
  i = 0
  with gz.open(outfile, 'wt') as out:
    with gz.open(infile, 'rt') as f:
      for line in f:
        if i == 0:
          spltln = line.split()
          out.write('\t'.join(spltln[0:4] + ['ID']) + '\n')
        else:
          spltln = line.split()
          if np.float32(spltln[4]) > 0.0:
            freq_values = np.array(spltln[4:], dtype=np.float32)
            # Generating the code
            test_values = simple_geodist_binning(freq_values, bins=bins)
            code = ''.join([str(i) for i in test_values])
            out.write('\t'.join(spltln[0:4] + [code]) + '\n')
        i += 1 

@click.command()
@click.option('--freqs', required=True, help='Minor allele count file')
@click.option('--bins', required=False, 
              type=str, default="[0.0, 0.01, 0.05]", 
              help='Bin endpoints for allele frequency')
def main(freqs, bins):
    #1. Generating the bins
    bin_endpts = list(map(float, bins.strip('[]').split(',')))
    bins = generate_bins(bin_endpts)
    eprint("Generated Allele Frequency Bins...")
    eprint(bins)
    
    #2. Reading through the frequency file line by line and printing it out
    with gz.open(freqs,'rt') as f:
      i = 0
      for line in f:
        if i == 0:
          spltln = line.split()
          print('\t'.join(spltln[0:5] + ['ID']))
        else:
          # 2a. parse the lines with everything before the 5th column as part of the file
          spltln = line.split()
          if np.float32(spltln[4]) > 0.0:
            freq_values = np.array(spltln[5:], dtype=np.float32)
            # Generating the code
            test_values = simple_geodist_binning(freq_values, bins=bins)
            code = ''.join([str(i) for i in test_values])
            print('\t'.join(spltln[0:5] + [code]))
        
        i += 1
        if i % 100000 == 0:
          eprint("Finished %d00 k Variants" % int(i/100000))
    
if __name__ =='__main__':
  main()
