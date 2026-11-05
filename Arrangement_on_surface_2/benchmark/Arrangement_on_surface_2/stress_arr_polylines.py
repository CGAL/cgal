#!/usr/bin/env python3
"""Stress tests the arrangement of polylines with `bench_arr_polylines --verify`.

Every case runs in its own process, so that crashes, assertion violations, and
hangs are caught like failed checks. The input of a failed case is saved, and
then reduced while the case keeps failing in the same way.

Build the benchmark against the tree under test, with the assertions enabled, e.g.:

  cmake -S . -B <build> -DCGAL_DIR=<cgal tree> -DCMAKE_BUILD_TYPE=Debug \\
    -DCMAKE_CXX_FLAGS="-O1 -fsanitize=address,undefined -fno-sanitize=vptr,bool -fno-omit-frame-pointer"
  cmake --build <build> --target bench_arr_polylines_rational

The vptr and bool checks are off; they reject the downcasts of the arrangement
iterators and the first call of Trapezoidal_decomposition_2::set_with_guarantees(),
and a rejection would end the case before the checks that follow.
Worthwhile configurations are the sanitizers with the rational kernel (above),
Epeck with -DCGAL_CHECK_EXPENSIVE, and -DCGAL_ALWAYS_LEFT_TO_RIGHT.

The exit status is 1 if a case failed.
"""

import argparse
import collections
import concurrent.futures
import os
import re
import subprocess
import sys
from pathlib import Path

# Tiny cases on coarse grids have dense overlaps and reduce well.
PRESETS = {
    "grid-walk": "-t walk -n 8 -p 12 --size 8 --step 2 --snap 1",
    "grid-monotone": "-t monotone -n 8 -p 12 --size 8 --step 2 --snap 1",
    "grid-random": "-t random -n 6 -p 8 --size 4 --step 1 --snap 1",
    "tiny-walk": "-t walk -n 3 -p 6 --size 3 --step 1.5 --snap 1",
    "dense-walk": "-t walk -n 30 -p 30 --size 15 --step 2 --snap 1",
    "dense-monotone": "-t monotone -n 30 -p 30 --size 15 --step 2 --snap 1",
    "fine-walk": "-t walk -n 20 -p 20 --size 100 --step 20",
    "fine-random": "-t random -n 10 -p 10 --size 100",
}

ENV = dict(os.environ,
           ASAN_OPTIONS="detect_leaks=0",
           UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1")


def run(cmd, timeout):
    """Runs the benchmark; returns a key of the failure or None if it passed."""
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=ENV)
    except subprocess.TimeoutExpired as e:
        checks = re.findall(r"^Verify (\S+):", (e.stdout or b"").decode(), re.MULTILINE)
        return "timeout after " + (checks[-1] if checks else "start")
    if res.returncode == 0:
        return None
    failed = re.search(r"^Verify (\S+): FAIL", res.stdout, re.MULTILINE)
    if failed:
        return "verify " + failed.group(1)
    sanitizer = re.search(r"^SUMMARY: (.*)$", res.stderr, re.MULTILINE)
    if sanitizer:
        return sanitizer.group(1)
    ubsan = re.search(r"^(\S+?):(\d+):\d+: runtime error: (.*)$", res.stderr, re.MULTILINE)
    if ubsan:
        error = re.sub(r"0x[0-9a-f]+|\b\d+\b", "<n>", ubsan.group(3))
        return f"{Path(ubsan.group(1)).name}:{ubsan.group(2)}: {error}"
    cgal = re.search(r"^File\s*:\s*(\S+)\s*^Line\s*:\s*(\d+)", res.stderr, re.MULTILINE)
    if cgal:
        return f"CGAL {Path(cgal.group(1)).name}:{cgal.group(2)}"
    return f"exit status {res.returncode}"


def read_input(filename):
    tokens = Path(filename).read_text().split()
    polylines, pos = [], 1
    for _ in range(int(tokens[0])):
        m = int(tokens[pos])
        coords = tokens[pos + 1 : pos + 1 + 2 * m]
        polylines.append(list(zip(coords[::2], coords[1::2])))
        pos += 1 + 2 * m
    return polylines


def write_input(filename, polylines):
    lines = [str(len(polylines))]
    for pl in polylines:
        lines.append(" ".join([str(len(pl))] + [c for pt in pl for c in pt]))
    Path(filename).write_text("\n".join(lines) + "\n")


def drop_repeats(points):
    """Removes consecutive repeats, which removing a point in between may create."""
    return [p for i, p in enumerate(points) if (i == 0) or (p != points[i - 1])]


def reduce_list(items, still_fails, min_size):
    """Removes chunks of items of decreasing sizes as long as the case fails."""
    chunk = max(len(items) // 2, 1)
    while True:
        pos, removed = 0, False
        while pos < len(items) and len(items) - chunk >= min_size:
            candidate = items[:pos] + items[pos + chunk:]
            if still_fails(candidate):
                items, removed = candidate, True
            else:
                pos += chunk
        if chunk == 1 and not removed:
            return items
        chunk = max(chunk // 2, 1)


def reduce_input(filename, reduced, base_cmd, key, timeout):
    trial = str(reduced) + ".trial"
    cmd = base_cmd + ["-i", trial]

    def fails(polylines):
        write_input(trial, polylines)
        return run(cmd, timeout) == key

    polylines = reduce_list(read_input(filename), fails, 1)
    for i in range(len(polylines)):
        def fails_with_points(points):
            points = drop_repeats(points)
            return len(points) >= 2 and fails(polylines[:i] + [points] + polylines[i + 1:])
        polylines[i] = drop_repeats(reduce_list(polylines[i], fails_with_points, 2))
    write_input(reduced, polylines)
    Path(trial).unlink(missing_ok=True)
    return polylines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bench", help="the bench_arr_polylines executable")
    parser.add_argument("-n", "--seeds", type=int, default=100, help="number of seeds per preset")
    parser.add_argument("--first-seed", type=int, default=0, help="the first seed")
    parser.add_argument("-p", "--preset", action="append", choices=sorted(PRESETS),
                        help="preset to run (can be repeated); by default, all")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="number of parallel cases")
    parser.add_argument("-t", "--timeout", type=float, default=300, help="time limit of a case in seconds")
    parser.add_argument("-o", "--output", default="stress_failures", help="directory of the failed inputs")
    parser.add_argument("--no-reduce", action="store_true", help="do not reduce the failed inputs")
    args = parser.parse_args()

    bench = str(Path(args.bench).resolve())
    presets = args.preset or list(PRESETS)
    cases = [(p, s) for p in presets for s in range(args.first_seed, args.first_seed + args.seeds)]

    def run_case(case):
        preset, seed = case
        return case, run([bench, "--verify", "-s", str(seed)] + PRESETS[preset].split(), args.timeout)

    failures = {}  # the first failed case of every key
    counts = collections.Counter()
    num_failed = 0
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as executor:
        for i, (case, key) in enumerate(executor.map(run_case, cases), 1):
            if key is not None:
                num_failed += 1
                print(f"FAILED {case[0]} seed {case[1]}: {key}", flush=True)
                failures.setdefault(key, case)
                counts[key] += 1
            if (i % 100 == 0) or (i == len(cases)):
                print(f"{i}/{len(cases)} cases, {num_failed} failed", flush=True)

    if not failures:
        return 0
    print()
    for key, count in counts.most_common():
        print(f"{count:6} {key}")

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    for key, (preset, seed) in failures.items():
        generator = [bench, "-s", str(seed)] + PRESETS[preset].split()
        input_file = output / f"{preset}-{seed}.txt"
        subprocess.run(generator + ["-m", "none", "-o", str(input_file)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ENV, check=True)
        print(f"\n{key}\n  input: {input_file}")
        if args.no_reduce:
            continue
        snap = re.search(r"--snap \S+", PRESETS[preset])
        replay = [bench, "--verify", "-s", str(seed)] + (snap.group(0).split() if snap else [])
        reduced_file = output / f"{preset}-{seed}.reduced.txt"
        if run(replay + ["-i", str(input_file)], args.timeout) != key:
            print("  does not reproduce from the input file; not reduced")
            continue
        polylines = reduce_input(input_file, reduced_file, replay, key, args.timeout)
        print(f"  reduced to {len(polylines)} polylines, {sum(map(len, polylines))} points: {reduced_file}")
        print("  replay: " + " ".join(replay + ["-i", str(reduced_file)]))
    return 1


if __name__ == "__main__":
    sys.exit(main())
