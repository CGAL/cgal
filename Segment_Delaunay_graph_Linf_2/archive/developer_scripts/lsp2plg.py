#!/usr/bin/env python

import sys

for line in sys.stdin:
    linesplit = line.strip().split()
    if len(linesplit) == 3:
        assert(linesplit[0] == 'p')
        print('1')
        print(linesplit[1] + ' ' + linesplit[2])
    elif len(linesplit) == 5:
        assert(linesplit[0] == 's')
        print('2')
        print(linesplit[1] + ' ' + linesplit[2])
        print(linesplit[3] + ' ' + linesplit[4])
    elif len(linesplit) == 0:
        print
