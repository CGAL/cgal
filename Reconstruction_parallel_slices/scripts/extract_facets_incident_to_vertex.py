#!/usr/bin/python

from sys import argv
from sys import exit
from sys import stderr

 


if len(argv) != 3:
  stderr.write("Usage: "+argv[0]+"file.off vertex_index\n")
  exit()

f=file(argv[1])

line=f.readline().split()
nbvertex=int(line[1])
nbfacet=int(line[2])

res=""
for i in range(0,nbvertex):
  res+=f.readline()

nbf=0
for i in range (0,nbfacet):
  line=f.readline()
  assert(len(line.split())==4)
  try:
    line.split().index(argv[2])
    res+=line
    nbf+=1
  except:
    ff=0
print "OFF ",nbvertex,nbf,0
print res
