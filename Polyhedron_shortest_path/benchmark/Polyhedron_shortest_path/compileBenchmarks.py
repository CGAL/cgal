#!/usr/bin/python

"""
The purpose of this file is to recompute the benchmarks section of the user manual.

Its intended usage is

python compileBenchmarks.py <tableFile> <randomseed>

Where:
- <tableFile> - the file to output the tables to
- <randomseed> - a positive integer to seed the randomizer

Upon completion, the contents of <tableFile> can be pasted 
directly into the user manual source file, containing the 
run benchmarks.

"""

import sys;
import glob;
import subprocess;
import tempfile;
import random;
import re;


class TestConfig:
  """
  Describes a benchmark test case to run over a set of models
  """
  
  def __init__(self, testName, testModels, numTrials, numSources, numQueries, rand):
    """
    testName - the name of the test as it should appear in the documentation file
    testModels - a list of .off model files to run the test on
    numTrials - the number of times the full test should be repeated on each mesh
    numSources - the number of source points to use in each trial
    numQueries - the number of query points to use in each trial
    rand - a python random number generator, used to pass a random key to the tests
    """
    self.testName = testName;
    self.testModels = testModels;
    self.numTrials = numTrials;
    self.numSources = numSources;
    self.numQueries = numQueries;
    self.rand = rand;
    
  def __str__(self):
    return "%s : #Trials = %d, #Sources = %d, #Queries = %d" % (self.testName, self.numTrials, self.numSources, self.numQueries);

def prepare_program():
  cmakeCall = subprocess.call(['cmake', '-DCMAKE_CXX_FLAGS=-O2 -Wall', '-DCMAKE_BUILD_TYPE=Debug', '.']);
  if cmakeCall != 0:
    return False;
  makeCall = subprocess.call(['make']);
  if makeCall != 0:
    return False;
  return True;
    
def run_benchmarks(testConfig, outputFile):
  print("Config: " + str(testConfig));
  for model in testConfig.testModels:
    sys.stdout.write("Model = %s ... " % model);
    sys.stdout.flush();
    result = subprocess.call(['./benchmark_shortest_paths.exe', '-p', model.strip(), '-r', str(testConfig.rand.randrange(65536)), '-t', str(testConfig.numTrials), '-n', str(testConfig.numSources), '-q', str(testConfig.numQueries)], stdout=outputFile);
    if result == 0:
      sys.stdout.write("Done.\n");
    else:
      sys.stdout.write("Error.\n");
  return True;

def process_file(file, infoset):
  currentFile = None;
  for line in file:
    if line.strip():
      tokens = list(map(lambda x: x.strip(), line.strip().split('|')));
      if len(tokens) >= 2:
        if tokens[0].lower() == "filename":
          currentFile = tokens[1];
          if currentFile not in infoset:
            infoset[currentFile] = {};
        else:
          infoset[currentFile][tokens[0].lower()] = tokens[1]
      else:
        print(line);

def print_table(infoSet, config, outFile):
  outFile.write("\subsection Polyhedron_shortest_pathBenchmark%s %s\n" % (re.sub(r"\s+", "", config.testName.strip()), config.testName));
  outFile.write("<center>\n");
  outFile.write("Model | Number of Vertices | Average Construction Time (s) | Average Query Time (s) | Peak Memory Usage (MB)\n");
  outFile.write("---|---|---|---|---\n");
  
  for key in sorted(infoSet.keys()):
    outFile.write(key.split('/')[-1] + " | " + 
      infoSet[key].get("num vertices", "<crashed>") + " | " +
      infoSet[key].get("construction", "<crashed>") + " | " +
      infoSet[key].get("query", "<crashed>") + " | " + 
      infoSet[key].get("memory (peak)", "<crashed>") + "\n");
  outFile.write("</center>\n");
  outFile.write('\n');


  
# Specify a file to output to, and optionally a random seed to ensure consistent tests are run
if len(sys.argv) <= 2:
  print("Usage: python compileBenchmarks.py <inputModels> <outputFile> [randomSeed]" % sys.argv[0]);
  sys.exit(0);
  
testModelsFile = open(sys.argv[1], "r");
testModels = list(testModelsFile);

outFile = open(sys.argv[2], "w");
  
if len(sys.argv) > 3:
  globalRand = random.Random(int(sys.argv[3]));
else:
  globalRand = random.Random();

# Here I create two benchmark sets, one using a single source point each, and one using 10 source points
testConfigs = [
  TestConfig("Single Source Point", testModels, 20, 1, 100, globalRand),
  TestConfig("Ten Source Points", testModels, 20, 10, 100, globalRand), 
];

if not prepare_program():
  sys.exit(1);
  
for config in testConfigs:
  infoFile = tempfile.TemporaryFile();
  run_benchmarks(config, infoFile);
  infoFile.seek(0);
  infoSet = {};
  process_file(infoFile, infoSet);
  infoFile.close();
  print_table(infoSet, config, outFile);
