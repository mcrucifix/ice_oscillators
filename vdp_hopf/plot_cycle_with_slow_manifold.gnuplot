set term svg 
set output 'cycle_vdp.svg'
f(x) = x**3/3 - x
set xrange [-2.2:2.2]
set nokey
plot 'cycle_vdp.dat' using 2:3 with lines lw 3, 'vdp_fixed_point.dat' using 2:3 ps 4 pt 6, f(x) lw 3

