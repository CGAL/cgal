import sys
import random

if len(sys.argv) != 1:
    print 'Usage:', sys.argv[0]
    sys.exit(1)

previous = 0

saveline = ""

for line in sys.stdin:
    if previous == 1:
        print "s", saveline.rstrip(), line.rstrip()
        previous = 0
    else:
        r = random.randint(0,1)
        if r == 0:
            print "p", line.rstrip()
            previous = 0
        else:
            previous = 1
            saveline = line

if previous == 1:
    print "p", saveline.strip()
