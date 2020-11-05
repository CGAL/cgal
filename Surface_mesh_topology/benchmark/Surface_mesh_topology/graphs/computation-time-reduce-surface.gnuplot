# set terminal postscript eps color 20 lw 3
# set output '| epstopdf -f -o=computation-time-reduce-surface.pdf'

set terminal svg fname 'Verdana' lw 2 # size 640 480 fname 'Verdana' lw 3
set output 'computation-time-reduce-surface.svg'

set key autotitle columnheader

set ylabel "Time (sec)"
set xlabel "#Darts"
set key left

# set xtics (0, '' 1, 2, 4, 8, 16, 32, 64)
#set ytics (0, 30, 60, 90, 120, 150, 180,200)

set xrange [0:30500000]
set yrange [0:35]
# set auto y

# set logscale x
# set logscale y

# set xtics ('5,000,000' 5000000, '12,000,000' 12000000, '19,000,000' 19000000, '26,000,000' 26000000, )
# set xtics ('5,000,000' 5000000, '11,000,000' 11000000, '17,000,000' 17000000, '23,000,000' 23000000, '29,000,000' 29000000 )
set xtics ('5,000,000' 5000000, '10,000,000' 10000000, '15,000,000' 15000000, '20,000,000' 20000000, '25,000,000' 25000000,  '30,000,000' 30000000 )

# set auto x
set format x '%.0f'

FIT_LIMIT=1.e-14
# f(x) = a*x**2 + b*x + c
f(x) = a*x + b
fit f(x) 'computation-time-reduce-surface.dat' using 1:9 via a,b
# fit f(x) 'computation-time-reduce-surface.dat' using (log($1)):(log($9)) via a,b

plot 'computation-time-reduce-surface.dat' using 1:9 with points title "Reduce surface computation", f(x) title 'Model Fit'

