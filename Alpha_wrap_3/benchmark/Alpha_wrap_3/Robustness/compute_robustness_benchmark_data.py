# Copyright (c) 2019-2023 Google LLC (USA).
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org).
#
# $URL$
# $Id$
# SPDX-License-Identifier: GPL-3.0-or-later
#
#
# Author(s)     : Pierre Alliez
#                 Michael Hemmer
#                 Cedric Portaneri
#
#!/usr/bin/python

import os, sys, subprocess, datetime, time, signal, getopt

def signal_handler(signum, frame):
  raise Exception("Timed out!")

def compute_robustness_benchmark_data(execname, filename, alpha, max_time):

  exit_codes = {
   0 : "VALID_SOLID_OUTPUT",
   1 : "INPUT_IS_INVALID",
   2 : "OUTPUT_IS_NOT_TRIANGLE_MESH",
   3 : "OUTPUT_IS_COMBINATORIAL_NON_MANIFOLD",
   4 : "OUTPUT_HAS_BORDERS",
   5 : "OUTPUT_HAS_DEGENERATED_FACES",
   6 : "OUTPUT_HAS_GEOMETRIC_SELF_INTERSECTIONS",
   7 : "OUTPUT_DOES_NOT_BOUND_VOLUME",
   8 : "OUTPUT_DOES_NOT_CONTAIN_INPUT",
   9 : "OUTPUT_DISTANCE_IS_TOO_LARGE",
   10 : "SIGSEGV",
   11 : "SIGABRT",
   12 : "SIGFPE",
   13 : "TIMEOUT"
  }

  exit_code = 0
  output = ""
  cmd = ("/usr/bin/time", "-v",
         execname, "-i",
         filename, "-a", alpha)
  proc = subprocess.Popen(
      cmd,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      start_new_session=True)

  try:
    outs, errs = proc.communicate(timeout=int(max_time))
    exit_code = proc.returncode
    output = outs.decode("utf-8") + errs.decode("utf-8")

    for output_line in output.split("\n"):
      if output_line == "Command terminated by signal 11":
        exit_code = 10
        continue
      elif output_line == "Command terminated by signal 6":
        exit_code = 11
        continue
      elif output_line == "Command terminated by signal 8":
        exit_code = 12
        continue

  except subprocess.TimeoutExpired:
    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    exit_code = 13
    output = "process ran too long"

  print(exit_codes[exit_code])

def main(argv):
  execname=""
  filename=""
  alpha=""
  max_time=""
  try:
    opts, args = getopt.getopt(sys.argv[1:], 'e:i:a:t:')
  except getopt.GetoptError:
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-e":
      execname = arg
    elif opt == "-i":
      filename = arg
    elif opt == "-a":
      alpha = arg
    elif opt == "-t":
      max_time = arg

  compute_robustness_benchmark_data(execname, filename, alpha, max_time)

if __name__ == "__main__":
  main(sys.argv[1:])
