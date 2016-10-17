import argparse
import fileinput
import sys
import io

parser = argparse.ArgumentParser(
    description='Translate a plg file by the given x, y arguments.')
parser.add_argument('x',
                    help='x coordinate of translation')
parser.add_argument('y',
                    help='y coordinate of translation')
args = parser.parse_args()
x = int(args.x)
y = int(args.y)

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
        print(linesplit[0])
        newpolygon = False
    else:
        assert(len(linesplit) == 2)
        print(str(int(linesplit[0]) + x) + ' ' +
              str(int(linesplit[1]) + y)        )
        remaining = remaining - 1
        if remaining == 0:
            newpolygon = True

assert(newpolygon == True)
