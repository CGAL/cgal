from sys import argv
from sys import stderr
import codecs
import re

f=codecs.open(argv[1], encoding='utf-8')
res=""
np_name="np"
while True:
  line = f.readline()
  if not line:
    break;
  res+=line

  #possible np name
  if re.search("bgl_namedparameters", line):
    m = re.search("[@\\\]param\s+([^ ]+)\s", line)
    if m:
      np_name=m.group(1)
      print("found "+np_name)

  if re.search("cgalNamedParamsBegin", line):
    params=[]
    while True:
      line = f.readline()
      if not line:
        break
      res+=line
      if re.search("cgalNamedParamsEnd", line):
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
        stderr.write("Function with several nps! Not handled yet\n")
        exit(1)
      m = re.search("^(\s*).*{", line)
      if m:
        s=m.group(1)+"  CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS("+np_name
        for p in params:
          s+=", "+p+"_t"
        s+=");"

        add_macro=True

        #do nothing if the macro is already there
        while True:
          line = f.readline()
          if not line:
            stderr.write("Error parsing the function\n")
            exit(1)
          if re.search("CGAL_CHECK_AUTHORIZED_NAMED_PARAMETERS", line):
            add_macro=False
            res+=line
            break
          else:
            if line=="\n":
              res+=line
            else:
              break # do we want to support comments?

        if add_macro:
          res+=s+"\n"
        break;

o=open("/tmp/out.h", 'w')
o.write(res)
o.close()

