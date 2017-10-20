from __future__ import print_function
from collections import defaultdict
import xml.etree.ElementTree as ET
import os
import errno

xml = open("Test.xml", 'rb').read();

def open_file_create_dir(filename, mode):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    return open(filename, mode)

root=ET.fromstring(xml)
testing = root.find('Testing')
testlist = testing.find('TestList')
nb = 0
tests_ids = {}
for t in testlist:
    tests_ids[t.text] = nb
    nb += 1

tests = {}
for t in testing.findall('Test'):
    tests[tests_ids[t.find('FullName').text]] = \
         { \
           "Name":   t.find('Name').text, \
           "Status": t.attrib['Status'], \
           "Output": t.find('Results').find('Measurement').find('Value').text, \
           "Labels": [l.text for l in t.find('Labels').findall('Label')] if t.find('Labels') is not None else ['UNKNOWN_LABEL'], \
         }

tests_per_label = defaultdict(list)
for t_id in range(0, len(tests)):
    t = tests[t_id]
    for label in t['Labels']:
        tests_per_label[label].append(t)

for label, tests in tests_per_label.items():
    with open_file_create_dir("{}/error.txt".format(label), 'w') as error:
        for t in tests:
            print("   {} of {}".format("successful " if (t['Status'] == 'passed') else "ERROR:     ", t['Name']), file=error)
            with open("{}/ProgramOutput.{}".format(label, t['Name']), 'w') as f:
                f.write(t['Output'] if t['Output'] != None else "")
