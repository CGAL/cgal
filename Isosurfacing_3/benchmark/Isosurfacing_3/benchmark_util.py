import os
import sys
import re
import subprocess
import pandas as pd

def print_stream(stream):
    while True:
        line = stream.readline()
        if not line:
            break
        print(line, end="")

def run(cmd, output=True):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    exit_code = process.wait()
    if output:
        print_stream(process.stdout)
    if exit_code != 0:
        print_stream(process.stderr)
        sys.exit(exit_code)
    return process

def build(scenario, kernel, algorithm, tag):
    run(["cmake", "-E", "make_directory", "build"])
    run(["cmake", "-B", "build", "-DCMAKE_BUILD_TYPE=Release", "-DSCENARIO=" + scenario, "-DKERNEL=" + kernel, "-DALGO=" + algorithm, "-DTAG=" + tag, "-DCGAL_DIR=../../../"])
    run(["make", "-C", "build"])

def execute(n, threads, times=1):
    time = 0
    for i in range(times):
        process = run(["likwid", "-g", "MEM_DP", "-c", "S0:0-" + str(threads - 1), "./build/benchmark", "-N", str(n)], False)

        for line in process.stdout.readlines():
            print(line)

            m = re.search(r'internal timer:\s*(\d*)', line)
            if m is not None:
                time += int(m.group(1))
        
    return time / times