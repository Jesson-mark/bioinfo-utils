#!/usr/bin/env python
import sys

print('\nseq \t length')
for i in range(1,len(sys.argv)):
  seq = sys.argv[i]
  print(seq, len(seq))

if len(sys.argv) == 3:
  print('Equal:', len(sys.argv[1]) == len(sys.argv[2]))

