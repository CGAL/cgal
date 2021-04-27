
reset

set terminal png size 800,500
set output 'nearest_neighbor_benchmark.png'

set grid
set style data lines
#set logscale x

set title 'Time to find the neighbors of a point using a tree'
set xlabel "Number of points"
set ylabel "Time (us)"
set key autotitle columnhead

set datafile separator ","
plot for [col=2:10] 'nearest_neighbor_benchmark.csv' using 1:col
