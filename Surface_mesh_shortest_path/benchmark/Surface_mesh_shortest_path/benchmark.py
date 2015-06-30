import subprocess;
import sys;

class TestConfig:
  """
  Describes a benchmark test case to run over a set of models
  """
  
  def __init__(self, testName, testModels, numTrials, numSources, numQueries, randSeed, kernel="epick"):
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
    self.randSeed = randSeed;
    self.kernel = kernel;
    
  def __str__(self):
    return "%s : Kernel = %s #Trials = %d, #Sources = %d, #Queries = %d" % (self.testName, self.kernel, self.numTrials, self.numSources, self.numQueries);

def read_all_lines(file):
  f = open(file, "r");
  lines = list(map(lambda l: l.strip(), f.readlines()));
  f.close();
  return lines;
    
def read_test_config(line):
  toks = line.strip().split(',');
  return TestConfig(toks[0], read_all_lines(toks[1]), int(toks[2]), int(toks[3]), int(toks[4]), int(toks[5]));

def prepare_program():
  """
  Call cmake for the benchmark program in the current directory.  Specifies O2 to ensure optimized build, but Debug to allow NDEBUG sections to be included (for memory benchmarking for example)
  """
  cmakeCall = subprocess.call(['cmake', '-DCMAKE_CXX_FLAGS=-O3 -Wall -DCGAL_NDEBUG', '-DCMAKE_BUILD_TYPE=Debug', '.']);
  if cmakeCall != 0:
    return False;
  makeCall = subprocess.call(['make']);
  if makeCall != 0:
    return False;
  return True;
    
def run_benchmarks(testConfig, outputFile):
  """
  Executes the benchmarks program for the given configuration.
  """
  print("Config: " + str(testConfig));
  for model in testConfig.testModels:
    sys.stdout.write("Model = %s ... " % model);
    sys.stdout.flush();
    result = subprocess.call(['./benchmark_shortest_paths.exe', '-k', testConfig.kernel, '-p', model.strip(), '-r', str(testConfig.randSeed), '-t', str(testConfig.numTrials), '-n', str(testConfig.numSources), '-q', str(testConfig.numQueries)], stdout=outputFile);
    if result == 0:
      sys.stdout.write("Done.\n");
    else:
      sys.stdout.write("Error.\n");
  return True;

def process_file(file, infoset):
  """
  Processes the benchmark output into a dictionary
  """
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
