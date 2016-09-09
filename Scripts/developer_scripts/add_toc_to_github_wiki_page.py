from sys import argv
from sys import exit
import codecs
import re

# a probably incomplete version to generate an anchor from a section name
def get_anchor(s):
  s = s.replace("`","")
  s = s.replace("(","")
  s = s.replace(")","")
  s = s.replace(".","")
  s = s.replace("#","")
  s = s.replace(":","")
  s = s.replace(",","")
  s = s.replace(";","")
  s = s.replace("/","")
  s = s.replace("<","")
  s = s.replace(">","")
  s = s.replace("+","")
  s = s.replace("=","")
  s = s.replace("?","")
  s = s.replace("@","")
  s = s.lstrip(" ")
  s = s.rstrip("\n")
  s = s.rstrip(" ")
  s = s.replace(" ","-")
  s = s.lower()
  return "#"+s

# indices the nesting level (first level allowed is ##)
def get_level(s):
  m = re.search('^(#+)\s', s)
  if m:
    return len(m.group(1))
  else:
    return 0

def get_name(s):
  m = re.search('^#+\s+(.*)\s*$', s)
  if m:
    return m.group(1)
  else:
    return "ERROR: Section name extraction"

#generate the entry for one section
def get_toc_entry(s):
  name = get_name(s)
  level = get_level(s)-2
  anchor = get_anchor(s)

  if level<0:
    return "ERROR: h1 section are not allowed"

  res="* ["+name+"]("+anchor+")"
  for i in range(0,level):
    res="  "+res
  return res

#now the main
if len(argv) < 2:
  print("Nothing done, no input file provided")
  exit()

input = argv[1]

f = codecs.open(input, 'r', encoding='utf-8')

if not f:
  print("Cannot open "+input+"\n")
  exit()

#look for <!--TOC--> the begin of the file
line=f.readline()
if line.find("<!--TOC-->")==-1:
  exit()

#skip current TOC
line=f.readline()
while line and line.find("<!--TOC-->")==-1:
  line=f.readline()

if not line:
  exit()

buffer=""
TOC="<!--TOC-->\n\n# Table of Contents\n"

verbatim_mode=False # to ignore verbatim mode while looking for sections
TOC_empty=True
for line in f.readlines():
  buffer+=line
  if verbatim_mode:
    if line[:3]=="```":
      verbatim_mode=False
  else:
    if line[:3]=="```":
      verbatim_mode=True
    else:
      if line[0]=="#":
        TOC+=(get_toc_entry(line)+"\n")
        TOC_empty=False
TOC+="\n<!--TOC-->\n"

if not TOC_empty:
  f.close()
  f = codecs.open(input, 'w', encoding='utf-8')
  f.write(TOC)
  f.write(buffer)
