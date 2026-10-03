#!/usr/bin/env python3
# Copyright (c) 2026 Tel-Aviv University (Israel).
# SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
"""Compare the CSV outputs of arr_traits_benchmark for two builds.

Each side accepts one or more CSV files; the per-repetition times of a
benchmark are pooled across files, so several interleaved runs can be
combined to average out drift (see run_comparison.sh).

For every benchmark the script reports the ratio candidate/baseline of the
median and of the minimum time, together with a 95% bootstrap confidence
interval for the ratio of medians. A benchmark is classified as

  regression   the ratio of medians exceeds 1 + threshold and the whole
               confidence interval lies above 1;
  improvement  the symmetric condition below 1;
  same         otherwise (no significant difference).

Benchmarks whose checksums differ between the two sides did not perform the
same computation; they are reported as MISMATCH and must be investigated
before their timings are trusted.

Exit status: 0 if there is no regression and no mismatch, 1 otherwise.
"""

import argparse
import csv
import random
import statistics
import sys


def read(files):
    """Return {name: (checksum, [times])} pooled over the given files."""
    data = {}
    for path in files:
        with open(path, newline="") as f:
            rows = csv.DictReader(line for line in f if not line.startswith("#"))
            for row in rows:
                times = [float(t) for t in row["times_ms"].split(";") if t]
                name = row["name"]
                checksum = row["checksum"]
                if name in data:
                    old_cs, old_times = data[name]
                    if old_cs != checksum:
                        sys.exit(f"{path}: checksum of {name} differs from "
                                 f"earlier files on the same side")
                    old_times.extend(times)
                else:
                    data[name] = (checksum, times)
    return data


def bootstrap_ratio_ci(base, cand, resamples, rng):
    """95% bootstrap confidence interval of median(cand) / median(base)."""
    ratios = []
    for _ in range(resamples):
        b = statistics.median(rng.choices(base, k=len(base)))
        c = statistics.median(rng.choices(cand, k=len(cand)))
        if b > 0:
            ratios.append(c / b)
    ratios.sort()
    if not ratios:
        return float("nan"), float("nan")
    lo = ratios[int(0.025 * (len(ratios) - 1))]
    hi = ratios[int(0.975 * (len(ratios) - 1))]
    return lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", nargs="+", required=True,
                    help="CSV file(s) of the baseline build (e.g., before the change)")
    ap.add_argument("--candidate", nargs="+", required=True,
                    help="CSV file(s) of the candidate build (e.g., after the change)")
    ap.add_argument("--threshold", type=float, default=5.0,
                    help="relative slowdown, in percent, considered significant "
                         "(default 5)")
    ap.add_argument("--resamples", type=int, default=2000,
                    help="bootstrap resamples (default 2000)")
    ap.add_argument("--seed", type=int, default=1,
                    help="seed of the bootstrap (default 1)")
    args = ap.parse_args()

    base = read(args.baseline)
    cand = read(args.candidate)
    rng = random.Random(args.seed)
    thr = args.threshold / 100.0

    header = (f"{'benchmark':30}{'base med':>11}{'cand med':>11}"
              f"{'ratio':>8}{'min ratio':>10}{'95% CI':>17}  verdict")
    print(header)
    print("-" * len(header))

    regressions = mismatches = 0
    for name in list(base) + [n for n in cand if n not in base]:
        if name not in base or name not in cand:
            side = "candidate" if name in cand else "baseline"
            print(f"{name:30}{'':>57}  only in {side}")
            continue
        b_cs, b_t = base[name]
        c_cs, c_t = cand[name]
        b_med, c_med = statistics.median(b_t), statistics.median(c_t)
        ratio = c_med / b_med if b_med > 0 else float("nan")
        min_ratio = min(c_t) / min(b_t) if min(b_t) > 0 else float("nan")
        lo, hi = bootstrap_ratio_ci(b_t, c_t, args.resamples, rng)
        if b_cs != c_cs:
            verdict = "MISMATCH (checksums differ)"
            mismatches += 1
        elif ratio > 1 + thr and lo > 1:
            verdict = "REGRESSION"
            regressions += 1
        elif ratio < 1 - thr and hi < 1:
            verdict = "improvement"
        else:
            verdict = "same"
        print(f"{name:30}{b_med:11.3f}{c_med:11.3f}{ratio:8.3f}{min_ratio:10.3f}"
              f"   [{lo:5.3f}, {hi:5.3f}]  {verdict}")

    print()
    print(f"{regressions} regression(s), {mismatches} checksum mismatch(es) "
          f"(threshold {args.threshold:g}%)")
    return 1 if regressions or mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
