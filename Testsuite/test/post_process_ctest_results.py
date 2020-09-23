import sys
import io
import re
import os

report_file=sys.argv[1]
report_name=sys.argv[2]
global_report_name=sys.argv[3]
rx=re.compile('(.*Configuring (examples|demo|test)*( in )*(test\/|examples\/|demo\/)*)((?!done)\w+)')
rx_demo=re.compile('.*in demo\/')
rx_examples=re.compile('.*in examples\/')


#open the Installation report
#For each NAME, check if NAME is a directory. If not, create one, create a 
#text report, and write everything that is in the report until the next NAME
#in it. Then, add 'NAME r' in the global report. This should allow to get all 
#the NOTICE and other info explaining why the configuration is skiped.

name=""
is_writing=False
is_ignored=False
global_report=open(global_report_name, "a+")
with open(report_file, "rt") as test_report:
  for myline in test_report:
    m=rx.match(myline)

    if is_writing:
      if m:
        is_writing=False
        test_report.close()
        if is_ignored:
          print("{label} {result}".format(label=name, result='r'), file=global_report)
          is_ignored=False
      else:
        test_report.write(myline)
    if not is_writing:
      if m:
        name=m.group(0).replace(m.group(1), "")
        if rx_demo.match(myline):
          name="{str}_Demo".format(str=name)
        elif rx_examples.match(myline):
          name="{str}_Examples".format(str=name)
        elif name == "libCGAL":
          name="libCGAL_shared"
        elif name == "libCGAL_Core":
          name="libCGALCore_shared"
        elif name == "libCGAL_ImageIO":
          name="libCGALimageIO_shared"
        elif name == "libCGAL_Qt5":
          name="libCGALQt5_shared"
        if name=="incomplete":
          is_writing=False
          is_ignored=False
          continue
        else:
          if not os.path.isdir(name):
            is_ignored=True
            os.mkdir(name)
            test_report=open("{dir}/{file}".format(dir=name, file=report_name), "w+")
            print("""
{scm_branch}
"""       .format(scm_branch=open("{}/../../../../../.scm-branch".format(os.getcwd()), 'r').read()),file=test_report)
          else:
            is_ignored=False
            test_report=open("{dir}/{file}".format(dir=name, file=report_name), "a+")
            test_report.write(" --- CMake Results: --- \n\n")
          is_writing=True
if is_writing:
    is_writing=False
    test_report.close()
    if is_ignored:
      print("{label} {result}".format(label=name, result='r'), file=global_report)
      is_ignored=False
global_report.close()
