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

params_stack={}
np_names=[]
while True:
  line = f.readline()
  if not line:
    break;
  res+=line

  if bool(params_stack):
    if not is_in_comment_section(line):
      m = re.search("^(\s*).*{", line)
      if m:
        s=""
        for np_name, params in params_stack.items():
          if len(params)==0:
            stderr.write("Cannot parse documented named parameters ("+argv[1]+")\n")
            exit(1)
          s+=m.group(1)+"  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS("+np_name
          for p in params:
            s+=", "+p+"_t"
          s+=");\n"
        params_stack={}

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
            continue # in case there are more than one macro
          else:
            #~ if line=="\n":
              #~ res+=line
            #~ else:
              break # do we want to support comments?

        res+=s+"\n"
        if line != "\n":
          res+=line
        modified=True
        continue;
  else:
    if not is_in_comment_section(line):
      continue

  #possible np name (TODO: does not work if not on the same line)
  if re.search("bgl_namedparameters", line):
    m = re.search("[@\\\]param\s+([^ ]+)\s", line)
    if m:
      np_name=m.group(1)
      np_names.append(np_name)

  if re.search("cgalNamedParamsBegin", line):
    params=[]
    while True:
      line = f.readline()
      if not line:
        break
      res+=line
      if re.search("cgalNamedParamsEnd", line):
        is_in_comment_section(line)
        for np_name in np_names:
          params_stack[np_name]=params
        np_names=[]
        break
      m = re.search(r"cgalParamNBegin{\s*([^ ]+)\s*}", line)
      if m:
        params.append(m.group(1).strip())
    if not line:
      stderr.write("ERROR: don't remember why I put that ("+argv[1]+")\n")
      break




f.close()
if modified:
  f=codecs.open(argv[1], encoding='utf-8', mode='w')
  f.write(res)
  f.close()
  #o=open("/tmp/out.h", 'w')
  #o.write(res)
  #o.close()

