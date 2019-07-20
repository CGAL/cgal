#!/usr/bin/env python

import sys

if len(sys.argv) != 1:
    print 'Usage:', sys.argv[0]
    sys.exit(1)

newpolygon = True
polysize = 0
remaining = 0
for line in sys.stdin:
    linesplit = line.split()
    if len(linesplit) == 0:
        continue
    if newpolygon == True:
        assert(len(linesplit) == 1)
        assert(linesplit[0].isdigit())
        polysize = int(linesplit[0])
        assert(polysize > 0)
        remaining = polysize
        newpolygon = False
    else:
        assert(len(linesplit) == 2)

        if remaining == polysize:
            # save first point of polygon
            savefirst = linesplit
            #print ('savefirst= ' + \
            #       savefirst[0] + ' ' + savefirst[1] )
        else:
            print ('s ' + previous[0]  + ' ' + previous[1] + ' ' \
                        + linesplit[0] + ' ' + linesplit[1] )

        previous = linesplit
        remaining = remaining - 1

        if remaining == 0:
            if polysize > 2:
                print ('s ' + linesplit[0] + ' ' + linesplit[1] + ' ' \
                            + savefirst[0] + ' ' + savefirst[1] )
            else:
                if polysize == 1:
                    print ('p ' + linesplit[0] + ' ' + linesplit[1])

            newpolygon = True

assert(newpolygon == True)
