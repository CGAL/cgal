#!/usr/bin/python

"""
The purpose of this file is to recompute the benchmarks section of the user manual.

Its intended usage is

python compileBenchmarks.py <modelsFile> <dataFileBase> <tableFileBase> <figureFileBase> <randomseed>

Where:
- <modelsFile> - a file, containing all models to run tests on
- <dataFileBase> - output data files will be of the form "<dataFileBase>_<modelname>.dat", where <modelname> is the basename of the model
- <tableFileBase> - output tables will be of the form "<outputTableBase>_#.txt", where # is the number of source points
- <figureFileBase> - output figures will be of the form "<outputBaseName>_<type>.png" where <type> is one of {query, construction, memory}
- <randomseed> - a positive integer to seed the randomizer

Upon completion, a set of data, table, and figure files will be generated over the testing models set.

"""

import sys;
import subprocess;
import tempfile;
import random;
import re;
import benchmark;
import os;

def make_table_file_name(tableFileBase, numSources):
  return tableFileBase + '_' + str(numSources) + ".txt";

def make_model_data_file_name(dataFileBase, modelFile):
  return dataFileBase + '_' + os.path.basename(modelFile).partition('.')[0] + ".dat";

def make_figure_file_name(figureFileBase, modelFile):
  return figureFileBase + '_' + os.path.basename(modelFile).partition('.')[0] + ".png";
  
def print_to_datafiles(config, dataFileName, modelInfo):
  file = open(dataFileName, "a");
  file.write("%d %d %s %s %s\n" % (config.numSources, modelInfo.get("num vertices", 0), modelInfo.get('construction', '0'), modelInfo.get('query', '0'), modelInfo.get('memory (peak)', '0')));
  file.close();

def print_table(infoSet, config, outFile):
  outFile.write("\subsection Polyhedron_shortest_pathBenchmark%s %s\n" % (re.sub(r"\s+", "", config.testName.strip()), config.testName));
  outFile.write("<center>\n");
  outFile.write("Model | Number of Vertices | Average Construction Time (s) | Average Query Time (s) | Peak Memory Usage (MB)\n");
  outFile.write("---|---|---|---|---\n");

  for key, value in sorted(infoSet.items(), key=lambda v: int(v[1]["num vertices"])):
    outFile.write(key.split('/')[-1] + " | " + 
      value.get("num vertices", "(crashed)") + " | " +
      value.get("construction", "(crashed)") + " | " +
      value.get("query", "(crashed)") + " | " + 
      value.get("memory (peak)", "(crashed)") + "\n");
  outFile.write("</center>\n");
  outFile.write('\n');
  
# Specify a file to output to, and optionally a random seed to ensure consistent tests are run
if len(sys.argv) <= 5:
  print("Usage: python %s <modelsFile> <dataFileBase> <tableFileBase> <figureFileBase> <randomseed>" % sys.argv[0]);
  sys.exit(0);
  
sampleRange = reversed([1] + list(range(5, 55, 5)));
  
testModels = benchmark.read_all_lines(sys.argv[1]);

if not benchmark.prepare_program():
  sys.exit(1);
  
rand = random.Random(int(sys.argv[5]));
  
dataFileBase = sys.argv[2];
tableFileBase = sys.argv[3];
figureFileBase = sys.argv[4];

for model in testModels:
  modelFileName = make_model_data_file_name(dataFileBase, model);
  if os.path.exists(modelFileName):
    os.remove(modelFileName);

for numSources in sampleRange:
  testname = "1 Source Point" if numSources == 1 else ("%d Source Points" % numSources);
  config = benchmark.TestConfig(testname, testModels, 20, numSources, 100, rand.randint(0, 65536));
  infoSet = {};
  infoFile = tempfile.TemporaryFile();
  benchmark.run_benchmarks(config, infoFile);
  infoFile.seek(0);
  infoSet = {};
  benchmark.process_file(infoFile, infoSet);
  infoFile.close();
  for k in infoSet.keys():
    print_to_datafiles(config, make_model_data_file_name(dataFileBase, k), infoSet[k]);
  tableFile = open(make_table_file_name(tableFileBase, numSources), "w");
  print_table(infoSet, config, tableFile);
  tableFile.close();
    
for runParams in [('query', 4, "Average Query Time"), ('construction', 3, "Average Construction Time"), ('memory', 5, "Peak Memory Usage")]:
  plotCommands = [];
    
  for k in testModels:
    plotCommands.append('"%s" using 1:%d with lines title "%s"' % (make_model_data_file_name(dataFileBase, k), runParams[1], os.path.basename(k)));

  plotCommand = "plot " + ', '.join(plotCommands);

  plotCommandFile = tempfile.TemporaryFile();
  plotCommandFile.write('set terminal png size 1280,960;\n');
  plotCommandFile.write('set output "%s";\n' % make_figure_file_name(figureFileBase, runParams[0]));
  plotCommandFile.write('set xlabel "Number of Source Points";\n');
  plotCommandFile.write('set ylabel "%s";\n' % runParams[2]);
  plotCommandFile.write(plotCommand + ";\n");
  plotCommandFile.write('unset output;\n');
  plotCommandFile.write('quit\n');
  plotCommandFile.seek(0);

  subprocess.call(['gnuplot'], stdin=plotCommandFile);
