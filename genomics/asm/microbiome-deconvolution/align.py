#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
from os import listdir
from os.path import isfile, join

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--folders', required=True)
parser.add_argument('-r', '--ref', required=True)
parser.add_argument('-a', '--align-dir', required=True)
parser.add_argument('-b', '--bwa-path', required=True)

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
  out_dir = args.align_dir + '/' + out_name
  files_to_merge = []

  cmd = ""
  cmd_prefix = '%s mem -t 8 %s' % (args.bwa_path, args.ref)

  for i, (f1, f2) in enumerate(pe_libs):
    cmd += '(' + cmd_prefix
    cmd += ' %s %s' % (folder+'/'+f1, folder+'/'+f2)
    cmd += ' | samtools view -bS - | samtools sort - %s.p%d ) &&' % (out_name, i)
    files_to_merge.append('%s.p%d.bam' % (out_name, i))
  for i, f in enumerate(se_libs):
    cmd += '( ' + cmd_prefix
    cmd +=' %s' % folder+'/'+f
    cmd += ' | samtools view -bS - | samtools sort - %s.u%d ) &&' % (out_name, i)
    files_to_merge.append('%s.u%d.bam' % (out_name, i))

  cmd += ' samtools merge %s.bam %s.*.bam' % (out_name, out_name)

  print cmd


