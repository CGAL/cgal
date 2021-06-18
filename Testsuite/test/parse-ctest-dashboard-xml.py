from __future__ import print_function
from collections import defaultdict
import xml.etree.ElementTree as ET
import os
import io
import errno
import re
import sys
import base64
import zlib
import gzip

result_file_name='{dir}/results_{tester}_{platform}.txt'
result_info_file_name='{dir}/results_{tester}_{platform}.info'
test_report_filename='{dir}/TestReport_{tester}_{platform}'

xml = open("Test.xml", 'rb').read()

def open_file_create_dir(filename, mode_, *args, **kwargs):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    if kwargs.get('gzip', None) == True:
        return gzip.open(filename, mode=mode_, encoding="utf-8")
    else:
        return io.open(filename, mode=mode_, encoding="utf-8")

root=ET.fromstring(xml)
testing = root.find('Testing')
testlist = testing.find('TestList')
nb = 0
tests_ids = {}
for t in testlist:
    tests_ids[t.text] = nb
    nb += 1

tests = {}
labels = set()
tester_name=sys.argv[1]
platform_name=sys.argv[2]
for t in testing.findall('Test'):
    nm = t.find('Results')
    t_exit_value=""
    ex_time = ""
    for ec in nm.findall("*/[@name='Execution Time']"):
      ex_time=ec.find('Value').text
    for ec in nm.findall("*/[@name='Exit Code']"):
      t_exit_value=ec.find('Value').text
    #t_exit_value=nm.find("*/[@name='Exit Code']").find('Value').text
    t_output = nm.find('Measurement').find('Value')
    t_output_value = t_output.text
    if t_output_value != None:
        if 'encoding' in t_output.attrib and t_output.attrib['encoding'] == 'base64':
           t_output_value = base64.standard_b64decode(t_output_value)
        if 'compression' in t_output.attrib and t_output.attrib['compression'] == 'gzip':
           t_output_value = zlib.decompress(t_output_value).decode("utf-8")
    tests[tests_ids[t.find('FullName').text]] = \
         { \
           "Name":   t.find('Name').text, \
           "Status": t.attrib['Status'], \
           "Output": t_output_value, \
           "Labels": [l.text for l in t.find('Labels').findall('Label')] if t.find('Labels') is not None else ['UNKNOWN_LABEL'], \
           "ExitValue": t_exit_value, \
           "ExecutionTime": ex_time, \
         }

tests_per_label = defaultdict(list)
for t_id in range(0, len(tests)):
    t = tests[t_id]
    for l in t['Labels']:
        if "_Tests" in l or "_Examples" in l or "_Demo" in l:
            label = l.replace("_Tests","")
            labels.add(label)
            tests_per_label[label].append(t)

warning_pattern=re.compile(r'(.*([^a-zA-Z_,:-])([^\d]\s)warning).*?(\[|\n)', flags=re.IGNORECASE)
w_det=re.compile("warning");
filter_pattern=re.compile(r'cmake|cgal', flags=re.IGNORECASE);
with open_file_create_dir(result_file_name.format(dir=os.getcwd(),
                                                  tester=tester_name,
                                                  platform=platform_name), 'a+') as results:
    for label, tests in tests_per_label.items():
        result_for_label='y'
        with open_file_create_dir("{}/error.txt".format(label), 'w') as error:
            for t in tests:
                print("   {result} {name} in {time} s : {value} ".format(result = "successful " if (t['Status'] == 'passed') else "ERROR:     ", name = t['Name'], value = t['ExitValue'] if(t['ExitValue'] != "") else "SUCCESS" , time = t['ExecutionTime']), file=error)
                if t['Status'] != 'passed':
                    result_for_label='n'
                elif t['Output'] != None and w_det.search(t['Output']):
                    entries = re.split("\n+", t['Output'])
                    for entry in entries:
                        m=warning_pattern.search(entry)
                        if m:
                            n = filter_pattern.search(m.group(0))
                            if n:
                                result_for_label='w'
                                break;
                            else:
                                result_for_label='t'

                with io.open("{}/ProgramOutput.{}".format(label, t['Name']), mode="w", encoding="utf-8") as f:
                    print("{}/ProgramOutput.{}".format(label, t['Name']))
                    f.write(t['Output'] if t['Output'] != None else "")


            print("{label} {result}".format(label=label, result=result_for_label), file=results)

for label, tests in tests_per_label.items():
        with open_file_create_dir(test_report_filename.format(dir=label,
                                                              tester=tester_name,
                                                              platform=platform_name), 'w') as label_report:
            print("""
{scm_branch}
------------------------------------------------------------------
- Error output from platform {platform}
------------------------------------------------------------------
{error_txt}
"""               .format(scm_branch=open("{}/../../../../../.scm-branch".format(os.getcwd()), 'r').read(),
                         platform=platform_name,
                         error_txt=open("{}/error.txt".format(label), 'r').read()), file=label_report)
            for t in tests:
                filename="{}/ProgramOutput.{}".format(label, t['Name'])
                with io.open(filename, mode="r", encoding="utf-8") as f:
                    print("""
------------------------------------------------------------------
- {file}
------------------------------------------------------------------

{file_text}
"""               .format(file=filename,
                          file_text=f.read()), file=label_report)
