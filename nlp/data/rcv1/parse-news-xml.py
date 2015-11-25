#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import xml.etree.ElementTree as et

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--filename', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

# input = open(args.filename)
input = open(sys.stdin)
text = '<data>\n' + input.read() + '</data>\n'
root = et.fromstring(text)

def isokay(sentence):
  words = sentence.split()
  if words[0].isupper():
    return False
  elif not words[0][0].isalpha():
    return False
  else:
    return True

def clean_up(text):
  words = [w for w in text.split() if not w.isupper() and w != '-']
  if len(words) <= 4:
    return ''
  else:
    return ' '.join(words)

headline_f = open(args.filename + '.1', 'w')
text_f     = open(args.filename + '.2', 'w')

for child in root:
  pair = [None, None]
  # print child.tag, child.attrib
  for el in child:
    # print el.tag, el.attrib
    if el.tag == 'headline':
      if el.text.startswith('PRESS DIGEST'): break
      pair[0] = clean_up(el.text)
    if el.tag  == 'text':
      pars = [c for c in el if c.tag == 'p' and isokay(c.text)]
      if pars:
        text = pars[0].text
        if text[-1] in ('.', '?', '!'):
          pair[1] = text
          break
  if pair[0] and pair[1]:
    # print pair
    headline_f.write(pair[0] + '\n')
    text_f.write(pair[1] + '\n')
