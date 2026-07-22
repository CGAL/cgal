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
    measurements = {"time" : 0, "polygons" : 0, "points" : 0, "bandwidth" : 0, "transfer" : 0, "performance" : 0, "clock" : 0, "intensity" : 0}

    for i in range(times):
        process = run(["likwid-perfctr", "-g", "MEM_DP", "-C", "S0:0-" + str(threads - 1), "./build/benchmark", "-N", str(n)], False)

        for line in process.stdout.readlines():
            print(line, end='')

            m = re.search(r'internal timer:\s*(\d+)', line)
            if m is not None:
                measurements["time"] += int(m.group(1))

            m = re.search(r'internal polygons:\s*(\d+)', line)
            if m is not None:
                measurements["polygons"] += int(m.group(1))

            m = re.search(r'internal points:\s*(\d+)', line)
            if m is not None:
                measurements["points"] += int(m.group(1))

            m = re.search(r'Memory bandwidth.*\s+(\d+(\.\d+)?) \|\s*$', line)
            if m is not None:
                measurements["bandwidth"] += float(m.group(1))

            m = re.search(r'Memory data volume.*\s+(\d+(\.\d+)?) \|\s*$', line)
            if m is not None:
                measurements["transfer"] += float(m.group(1))

            m = re.search(r'DP.*\s+(\d+(\.\d+)?) \|\s*$', line)
            if m is not None:
                measurements["performance"] += float(m.group(1))

            m = re.search(r'Clock.*\s+(\d+(\.\d+)?) \|\s*$', line)
            if m is not None:
                measurements["clock"] += float(m.group(1))

            m = re.search(r'Operational intensity.*\s+(\d+(\.\d+)?) \|\s*$', line)
            if m is not None:
                measurements["intensity"] += float(m.group(1))

    for item in measurements.items():
        measurements[item[0]] = item[1] / times

    return measurements