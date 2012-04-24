#!/usr/bin/python
#this dirty scripts tries to put one copyright holder per line
#so that it is easy to get them out from a file
#WARNING: this is an experimental script that is not guarantee to work
# in all cases
#usage: cgal_header_clean_up.py cpp_file


from sys import argv
from sys import exit
from sys import stderr
import re
import curses.ascii

f=file(argv[1])

line=f.readline();


if re.search("// Copyright \(c\) ",line):
  index=16
  
  while (not curses.ascii.isalpha(line[index]) and line[index]!='\n' ):
    index=index+1
  print line[0:index]
  print "// "+line[index:len(line)-1]
 
  regexp = re.compile("//\s*$")

  line=f.readline()
  buf=""
  while 1:
    if regexp.match(line):
      break
    #look for a comma
    index=3
    last=3
    
    while 1:
      while (line[index]!=',' and line[index]!='\n'):
        index=index+1
      if line[index]=='\n':
        if last!=index:
          buf=line[last:index]+" "
        else:
          buf=""
        break
      print "// "+buf+line[last:index+1]
      buf=""
      last=index+1
      while line[last]==" ":
        last=last+1
      index=last
    line=f.readline()
  print "// "+buf
  while line:
    print line,
    line=f.readline()
  
    
else:
  stderr.write(argv[1]+" error no copyright in the first place\n")  
  exit()
