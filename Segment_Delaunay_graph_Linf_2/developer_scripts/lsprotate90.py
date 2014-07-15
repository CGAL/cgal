#!/usr/bin/env python

import sys

def inv(s):
  if s[0] == '-':
    return s[1:]
  elif s[0] == '+':
    return '-' + s[1:]
  else: # plain number
    return '-' + s

if len(sys.argv) != 1:
  print 'Usage:', sys.argv[0]
  sys.exit(1)

for line in sys.stdin:
  linesplit = line.strip().split()
  if len(linesplit) == 3:
    assert(linesplit[0] == 'p')
    print('p ' + inv(linesplit[2]) + ' ' + linesplit[1])
  elif len(linesplit) == 5:
    assert(linesplit[0] == 's')
    print('s ' + \
          inv(linesplit[2]) + ' ' + linesplit[1] + ' ' + \
          inv(linesplit[4]) + ' ' + linesplit[3] )
  elif len(linesplit) == 0:
    print
