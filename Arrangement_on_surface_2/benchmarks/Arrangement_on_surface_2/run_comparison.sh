#!/bin/sh
# Copyright (c) 2026 Tel-Aviv University (Israel).
# SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
#
# Run the baseline and candidate builds of arr_traits_benchmark alternately,
# so that slow drift of the machine (thermal throttling, background load)
# affects both builds equally, then compare the pooled results.
#
# Usage: run_comparison.sh <baseline-exe> <candidate-exe> [rounds] [-- bench options]
# Example:
#   ./run_comparison.sh build-old/arr_traits_benchmark \
#                       build-new/arr_traits_benchmark 5 -- --seed 7 -r 5
#
# For stable timings, run on an otherwise idle machine and, on Linux, pin the
# process to one core, e.g. by prefixing the command with "taskset -c 2".

set -e

if [ $# -lt 2 ]; then
  sed -n '6,15p' "$0"
  exit 2
fi

baseline=$1
candidate=$2
shift 2
rounds=3
if [ $# -gt 0 ] && [ "$1" != "--" ]; then
  rounds=$1
  shift
fi
[ "$1" = "--" ] && shift

dir=$(mktemp -d "${TMPDIR:-/tmp}/arr_bench.XXXXXX")
here=$(cd "$(dirname "$0")" && pwd)

i=1
while [ "$i" -le "$rounds" ]; do
  echo "round $i/$rounds: baseline" >&2
  "$baseline" "$@" --csv "$dir/baseline_$i.csv" > /dev/null
  echo "round $i/$rounds: candidate" >&2
  "$candidate" "$@" --csv "$dir/candidate_$i.csv" > /dev/null
  i=$((i + 1))
done

echo "results in $dir" >&2
python3 "$here/compare_benchmarks.py" \
  --baseline "$dir"/baseline_*.csv --candidate "$dir"/candidate_*.csv
