#!/usr/local/bin/python3

'''
  Naive categorization of geodist categories
'''

import numpy as np 
import gzip as gz

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
