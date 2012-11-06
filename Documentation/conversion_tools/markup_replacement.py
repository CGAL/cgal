#!/usr/bin/python2

#replace markup #, ## ,### by \section, \subsection, \subsubsection.
#anchor names are preserved and generated from the section name otherwise
#The script is not perfect and might miss some specific cases

from sys import argv
from os import path
import string
import re

anchors={}

def generate_anchor(chapter,text):
  pattern = re.compile('[\W_]+')
  words=text.split()
  i=1;
  res=chapter+pattern.sub('',words[0])
  while len(res)<40 and i<len(words):
    word=pattern.sub('',words[i])
    res+=word
    i+=1
  if anchors.has_key(res):
    anchors[res]+=1
    res+="_"+str(anchors[res])
  else:
    anchors[res]=0
  return res

f=file(argv[1])
regexp_line=re.compile('^\s*#')
#~ regexp_section=re.compile('^\s*#\s*([ a-b().,]+)\s*#(.*)')
regexp_section=re.compile('^\s*(#+)\s*([0-9a-zA-Z (),.:?%-`\']+[0-9a-zA-Z.?`)])\s*#+(.*)')
regexp_anchor=re.compile('^\s*{#([0-9a-zA-Z_]+)}')

result=""
diff=False
chapter=path.abspath(argv[1]).split('/')[-2]
for line in f.readlines():
  if regexp_line.match(line):
    m=regexp_section.search(line)
    if m:
      values=m.groups()
      anchor=''
      if len(values)==2:
        anchor=generate_anchor(chapter,values[1])
      else:
        anchor=regexp_anchor.match(values[2])
        if anchor:
          anchor=anchor.group(1)
        else:
          anchor=generate_anchor(chapter,values[1])
      if len(values[0])==1:
        result+="\section "+anchor+" "+values[1]+"\n"
      elif len(values[0])==2:
        result+="\subsection "+anchor+" "+values[1]+"\n"
      elif len(values[0])==3:
        result+="\subsubsection "+anchor+" "+values[1]+"\n"
      else:
        print "Error while processing "+argv[1]
        assert False
      diff=True
    else:
      result+=line
  else:
    result+=line
f.close()

if diff:
  f=file(argv[1],'w')
  f.write(result)
  f.close()
