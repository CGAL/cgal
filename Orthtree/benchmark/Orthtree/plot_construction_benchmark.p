
reset

set terminal png size 800,500
set output 'construction_benchmark.png'

set grid
set style data lines
#set logscale x

set title 'Time to construct a tree from points'
set xlabel "Number of points"
set ylabel "Time (ms)"
set key autotitle columnhead

set datafile separator ","
plot for [col=2:10] 'construction_benchmark.csv' using 1:col
