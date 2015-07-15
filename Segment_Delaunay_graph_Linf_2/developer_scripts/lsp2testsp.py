#!/usr/bin/env python

import sys

if len(sys.argv) != 1:
    print 'Usage:', sys.argv[0]
    sys.exit(1)

for line in sys.stdin:
    linesplit = line.strip().split()
    if len(linesplit) == 3:
        assert(linesplit[0] == 'p')
        print('Point_2(' + linesplit[1] + ', ' + linesplit[2] + '),')
    elif len(linesplit) == 5:
        assert(linesplit[0] == 's')
        print('Segment_2(' + \
              'Point_2(' + linesplit[1] + ', ' + linesplit[2] + '), ' + \
              'Point_2(' + linesplit[3] + ', ' + linesplit[4] + ')),' )
    elif len(linesplit) == 0:
        print
