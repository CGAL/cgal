set terminal postscript eps color
set output "benchmark.eps"
set xlabel "total number of boxes"
set ylabel "time spent in seconds"
set autoscale
set key left
set multiplot
set logscale xy
show ylabel
show xlabel

plot "benchmark.data"  title "Streaming" w l lw 2, \
      "benchmark.data" using 1:3 title "Scanning" w l lw 2


