#!/usr/bin/python2

#checks whether there is a duplicated anchor name used in \section, \subsection or \subsubsection

from sys import argv
import re

anchors={}

pattern=re.compile("\\\(sub)?(sub)?section\s+(\w*)")

for i in range(0,len(argv)):
  f=file(argv[i])
  for line in f.readlines():
    res=pattern.search(line)
    if res:
      key=res.groups()[-1]
      if anchors.has_key(key):
        print anchors[key],":",key, "in", argv[i]
      else:
        anchors[key]=argv[i]