# set terminal postscript eps color 20 lw 3
# set output '| epstopdf -f -o=computation-time-polygonal-schema.pdf'

set terminal svg fname 'Verdana' lw 2 # size 640 480 fname 'Verdana' lw 3
set output 'computation-time-polygonal-schema.svg'

set key autotitle columnheader

set ylabel "Time (sec)"
set xlabel "Path lengths"
set key left

# set xtics (0, '' 1, 2, 4, 8, 16, 32, 64)
# set ytics (4, 16, 64, 256, 1024, 4096, "16,384" 16384)

# set logscale x 10
# set logscale y 10

set xrange [0:62000000]
# set yrange [0:250]

# set xtics ('5,000,000' 5000000, '10,000,000' 10000000, '15,000,000' 15000000, '20,000,000' 20000000, '25,000,000' 25000000, '30,000,000' 30000000)
set xtics ('10,000,000' 10000000, '20,000,000' 20000000, '30,000,000' 30000000, '40,000,000' 40000000, '50,000,000' 50000000, '60,000,000' 60000000)

# set auto x
 
FIT_LIMIT=1.e-14
f(x) = m*x + b
fit f(x) 'computation-time-polygonal-schema.dat' using ($1+$2):4 via m,b

plot 'computation-time-polygonal-schema.dat' using ($1+$2):4 with points title "Homotopy test", f(x) title 'Model Fit'

