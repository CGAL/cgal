# Copyright (c) 2019-2023 Google LLC (USA).
# All rights reserved.
#
# This file is part of the 3D Alpha Wrapping package, which is being prepared for
# submission to CGAL (www.cgal.org).
#
# $URL$
# $Id$
# SPDX-License-Identifier: GPL-3.0-or-later
#
#
# Author(s)     : Pierre Alliez
#                 Michael Hemmer
#                 Cedric Portaneri
#                 Mael Rouxel-Labbé
#

#!/bin/bash

# $1: directory containing the alpha wrap project
# $2: directory containing the output results
# $3: alpha value
# $4: timeout value in seconds
# $5: hash of the latest commit
# $6: the input file path
function compute_benchmark_data() {
  filename=$(basename -- "$6")
  filename="${filename%.*}"

  python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Robustness/compute_robustness_benchmark_data.py \
   -e $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/build-release/robustness_benchmark -i $6 -a $3 -t $4 \
   > $2/Robustness/results/$5/$filename.log

  python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Performance/compute_performance_benchmark_data.py \
   -e $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/build-release/performance_benchmark -i $6 -a $3 \
   > $2/Performance/results/$5/$filename.log

  python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Quality/compute_quality_benchmark_data.py \
   -e $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/build-release/quality_benchmark -i $6 -a $3 \
   > $2/Quality/results/$5/$filename.log
}
export -f compute_benchmark_data

# $1: directory containing the alpha wrap project
# $2: directory containing the input data folder
# $3: directory containing the output results
# $4: alpha value
# $5: timeout value for robustness benchmark in seconds
# $6: number of virtual thread used
# $7: hash of the latest commit
# $8: hash of a commit to perform the diff with latest
cd $1

mkdir -p $3/Robustness/results/$7
mkdir -p $3/Performance/results/$7
mkdir -p $3/Quality/results/$7
mkdir -p $3/Robustness/charts_data
mkdir -p $3/Performance/charts_data
mkdir -p $3/Quality/charts_data
mkdir -p $3/Robustness/charts
mkdir -p $3/Performance/charts
mkdir -p $3/Quality/charts
mkdir -p $3/Robustness/log
mkdir -p $3/Performance/log
mkdir -p $3/Quality/log
mkdir -p $3/charts

find $2 -mindepth 1 | parallel -j$6 compute_benchmark_data $1 $3 $4 $5 $7 :::

python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Robustness/generate_robustness_benchmark_charts.py -i $3/Robustness/results/$7 -o $3/Robustness -a $4 -c $7

if [ -z "$8" ]; then
  python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Performance/generate_performance_benchmark_charts.py -i $3/Performance/results/$7 -o $3/Performance -a $4 -c $7;
else
  python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Performance/generate_performance_benchmark_charts.py -i $3/Performance/results/$7 -o $3/Performance -a $4 -c $7 -p $3/Performance/results/$8 -d $8;
fi

if [ -z "$8" ]; then
  python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Quality/generate_quality_benchmark_charts.py -i $3/Quality/results/$7 -o $3/Quality -a $4 -c $7;
else
  python3 $1/Alpha_wrap_3/benchmark/Alpha_wrap_3/Quality/generate_quality_benchmark_charts.py -i $3/Quality/results/$7 -o $3/Quality -a $4 -c $7 -p $3/Quality/results/$8 -d $8;
fi

charts_path="$(ls "$3/Robustness/charts"/* -dArt | tail -n 1) $(ls "$3/Performance/charts"/* -dArt | tail -n 2) $(ls "$3/Quality/charts"/* -dArt | tail -n 9)"

pdfjam --nup 2x6 $charts_path --outfile $3/charts/results_$7_$8_alpha_$4_$(date '+%Y-%m-%d_%H:%M:%S').pdf
