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
import argparse;

def make_table_file_name(tableFileBase, numSources):
  return tableFileBase + '_' + str(numSources) + ".txt";

def make_model_data_file_name(dataFileBase, modelFile):
  return dataFileBase + '_' + os.path.basename(modelFile).partition('.')[0] + ".dat";

def make_figure_file_name(figureFileBase, modelFile):
  return figureFileBase + '_' + os.path.basename(modelFile).partition('.')[0] + ".png";
  
def print_to_datafiles(config, dataFileName, modelInfo):
  file = open(dataFileName, "a");
  file.write("%d %d %s %s %s\n" % (config.numSources, int(modelInfo.get("num vertices", 0)), modelInfo.get('construction', '0'), modelInfo.get('query', '0'), modelInfo.get('memory (peak)', '0')));
  file.close();

def print_table(infoSet, config, outFile):
  outFile.write("\subsection Surface_mesh_shortest_pathBenchmark%s %s\n" % (re.sub(r"\s+", "", config.testName.strip()), config.testName));
  outFile.write("<center>\n");
  outFile.write("Model | Number of Vertices | Average Construction Time (s) | Average Queries Per Second | Peak Memory Usage (MB)\n");
  outFile.write("---|---|---|---|---\n");

  for key, value in sorted(infoSet.items(), key=lambda v: int(v[1]["num vertices"])):
    outFile.write(key.split('/')[-1] + " | " + 
      value.get("num vertices", "(crashed)") + " | " +
      value.get("construction", "(crashed)") + " | " +
      value.get("query", "(crashed)") + " | " + 
      value.get("memory (peak)", "(crashed)") + "\n");
  outFile.write("</center>\n");
  outFile.write('\n');
  
parser = argparse.ArgumentParser(description="Run benchmarks on multiple model files");

parser.add_argument('-r', '--range', type=str, default=['1'], nargs='*', help="Can specify multiple single values, or ranges. Format is '#' or '#,#' or '#,#,#' for a single value, a range of values, or a range values of values with a specified skip");
parser.add_argument('-f', '--modelsfile', '--modelsfiles', type=str, nargs='*', help="a file containing a list of models to test on, one per line");
parser.add_argument('-m', '--model', '--models', type=str, nargs='*');
parser.add_argument('-s', '--randseed', type=int, help="Random seed for tests.");
parser.add_argument('-d', '--datafilebase', type=str, default='_modeldata', help="Base name for temporary data files (existing files will be overwritten).");
parser.add_argument('-t', '--tablefilebase', type=str, help="Base name for generated user manual tables (will be of the form \"<tablefilebase>_#.txt\", leave blank to suppress generation)");
parser.add_argument('-o', '--plotfilebase', type=str, help="Base name for generated plot figures (will be of the form \"<plotfilebase>_<modelname>.png\", leave blank to suppress generation)");
parser.add_argument('-k', '--kernel', choices=['ipick', 'epick', 'epeck'], default='epick', help="Geometry kernel to use for the benchmark.");
parser.add_argument('-n', '--numtrials', type=int, default=20, help="Specify the number of complete trial runs on each model.");
parser.add_argument('-q', '--numqueries', type=int, default=100, help="Specify the number of shortest path queries to run per trial.");

programArgs = parser.parse_args();

setOfSamples = set();

if programArgs.range:
  for r in programArgs.range:
    toks = r.split(',');
    if len(toks) == 1:
      setOfSamples.add(int(toks[0]));
    elif len(toks) == 2:
      setOfSamples.update(range(int(toks[0]), int(toks[1])));
    elif len(toks) == 3:
      setOfSamples.update(range(int(toks[0]), int(toks[1]), int(toks[2])));
    else:
      raise "Invalid range: " + r;

sampleRange = reversed(sorted(list(setOfSamples)));

testModels = [];

if programArgs.modelsfile:
  for modelsFile in programArgs.modelsfile:
    testModels.extend(benchmark.read_all_lines(modelsFile));

if programArgs.model:
  for model in programArgs.model:
    testModels.append(model);
    
if len(testModels) == 0:
  print("Error, must specify at least one model to benchmark");
  sys.exit(1);

rand = random.Random(programArgs.randseed);

dataFileBase = programArgs.datafilebase;
tableFileBase = programArgs.tablefilebase;
plotFileBase = programArgs.plotfilebase;

kernel = programArgs.kernel;

if tableFileBase == None and plotFileBase == None:
  print("Error, must specify either a table output or figure output file");
  sys.exit(1);

if not benchmark.prepare_program():
  print("Error, could not compile program");
  sys.exit(1);

for model in testModels:
  modelDataFileName = make_model_data_file_name(dataFileBase, model);
  if os.path.exists(modelDataFileName):
    os.remove(modelDataFileName);

for numSources in sampleRange:
  testname = "1 Source Point" if numSources == 1 else ("%d Source Points" % numSources);
  config = benchmark.TestConfig(testname, testModels, programArgs.numtrials, numSources, programArgs.numqueries, rand.randint(0, 65536), kernel);
  infoSet = {};
  infoFile = tempfile.TemporaryFile();
  benchmark.run_benchmarks(config, infoFile);
  infoFile.seek(0);
  infoSet = {};
  benchmark.process_file(infoFile, infoSet);
  infoFile.close();
  for k in infoSet.keys():
    print_to_datafiles(config, make_model_data_file_name(dataFileBase, k), infoSet[k]);
  if tableFileBase != None:
    tableFile = open(make_table_file_name(tableFileBase, numSources), "w");
    print_table(infoSet, config, tableFile);
    tableFile.close();

if plotFileBase != None:
  for runParams in [('query', 4, "Average Queries Per Second"), ('construction', 3, "Average Construction Time"), ('memory', 5, "Peak Memory Usage")]:
    plotCommands = [];
    for k in testModels:
      plotCommands.append('"%s" using 1:%d with lines title "%s"' % (make_model_data_file_name(dataFileBase, k), runParams[1], os.path.basename(k)));
    plotCommand = "plot " + ', '.join(plotCommands);
    plotCommandFile = tempfile.TemporaryFile();
    plotCommandFile.write('set terminal png size 640,480;\n');
    plotCommandFile.write('set output "%s";\n' % make_figure_file_name(plotFileBase, runParams[0]));
    plotCommandFile.write('set xlabel "Number of Source Points";\n');
    plotCommandFile.write('set ylabel "%s";\n' % runParams[2]);
    plotCommandFile.write(plotCommand + ";\n");
    plotCommandFile.write('unset output;\n');
    plotCommandFile.write('quit\n');
    plotCommandFile.seek(0);
    subprocess.call(['gnuplot'], stdin=plotCommandFile);
