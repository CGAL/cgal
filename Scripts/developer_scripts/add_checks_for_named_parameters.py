from sys import argv
from sys import stderr
import codecs
import re

f=codecs.open(argv[1], encoding='utf-8')
res=""
np_name="np"

modified=False


def is_in_comment_section(line):
  global in_comment_section

  if re.search("^\s*"+re.escape("//"), line):
    return True
  if re.search("^\s*"+re.escape("/*"), line):
    in_comment_section=True
    return True
  if not in_comment_section:
    return False
  if re.search(re.escape("*/"), line):
    in_comment_section=False;
  return True

in_comment_section = False

while True:
  line = f.readline()
  if not line:
    break;
  res+=line

  if not is_in_comment_section(line):
    continue

  #possible np name (TODO: does not work if not on the same line)
  if re.search("bgl_namedparameters", line):
    m = re.search("[@\\\]param\s+([^ ]+)\s", line)
    if m:
      np_name=m.group(1)
      #print("found "+np_name)

  if re.search("cgalNamedParamsBegin", line):
    params=[]
    while True:
      line = f.readline()
      if not line:
        break
      res+=line
      if re.search("cgalNamedParamsEnd", line):
        is_in_comment_section(line)
        break
      m = re.search(r"cgalParamNBegin{\s*([^ ]+)\s*}", line)
      if m:
        params.append(m.group(1).strip())
    if not line:
      break
    while True:
      line = f.readline()
      if not line:
        break
      res+=line
      if re.search("cgalNamedParamsBegin", line):
        stderr.write("Function with several nps! Not handled yet ("+argv[1]+")\n")
        exit(1)
      if is_in_comment_section(line):
        continue

      m = re.search("^(\s*).*{", line)
      if m:
        s=m.group(1)+"  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS("+np_name
        for p in params:
          s+=", "+p+"_t"
        s+=");"

        mb = re.search("(^\s+{)(.*)}$", line)
        if mb:
          res+=mb.group(1)+"\n";
          res+=s+"\n"
          if mb.group(2).strip():
            res+=m.group(1)+"  "+mb.group(2)+"\n"
          res+=m.group(1)+"}\n"
          modified=True
          break


        add_macro=True
        macro_already_here=False

        #do nothing if the macro is already there
        while True:
          line = f.readline()
          if not line:
            stderr.write("Error parsing the function\n")
            exit(1)
          if re.search("CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS", line):
            macro_already_here=True
            break
          else:
            if line=="\n":
              res+=line
            else:
              break # do we want to support comments?

        if add_macro:
          res+=s+"\n"
          if not macro_already_here:
            res+="\n"+line
          modified=True
        break;

f.close()
if modified:
  f=codecs.open(argv[1], encoding='utf-8', mode='w')
  f.write(res)
  f.close()
  #o=open("/tmp/out.h", 'w')
  #o.write(res)
  #o.close()

