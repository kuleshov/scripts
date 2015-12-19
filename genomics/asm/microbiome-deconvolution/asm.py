#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
from os import listdir
from os.path import isfile, join

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--folders', required=True)
parser.add_argument('-a', '--asm-dir', required=True)
parser.add_argument('-s', '--spades-path', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

folders = []
with open(args.folders) as f:
  folders = [line.strip() for line in f]

for folder in folders:
  files = [ f for f in listdir(folder) if isfile(join(folder,f)) ]
  pe_libs, se_libs = [], []
  for f in files:
    if not f.endswith('.fastq.gz'):
      continue
    if f.endswith('_1.fastq.gz'):
      root_f = f[:-11]
      pe_libs.append((root_f + '_1.fastq.gz', root_f + '_2.fastq.gz'))
    elif f.endswith('_2.fastq.gz'):
      continue
    else:
      se_libs.append(f)

  out_name = folder.strip().strip('/').split('/')[-1]
  out_dir = args.asm_dir + '/' + out_name
  cmd = '%s -o %s' % (args.spades_path, out_name)
  for i, (f1, f2) in enumerate(pe_libs):
    cmd += ' --pe%d-1 %s --pe%d-2 %s' % (i+1, folder+'/'+f1, i+1, folder+'/'+f2)
  for i, f in enumerate(se_libs):
    cmd += ' --s%d %s' % (i+1, folder+'/'+f)

  print cmd
  os.system(cmd)

