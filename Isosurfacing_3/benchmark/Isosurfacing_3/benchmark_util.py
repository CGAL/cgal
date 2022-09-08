import os
import sys
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

def build(scenario, kernel, algorithm):
    run(["cmake", "-E", "make_directory", "build"])
    run(["cmake", "-B", "build", "-DCMAKE_BUILD_TYPE=Release", "-DCGAL_DIR=../../../", "-D" + scenario + "= -D" + kernel + "= -D" + algorithm + "="])
    run(["cmake", "--build", "build"])

def execute(n, times=1):
    time = 0
    for i in range(times):
        run(["./build/benchmark", "-N", str(n)], False)
        time += int(process.stdout.readline())
    return time / times