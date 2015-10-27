#!/usr/bin/env python
# -*- coding: utf-8 -*-

from libkuleshov.stats import n50
import sys

L=[float(x) for x in sys.stdin]
print sum(L) / len(L)
