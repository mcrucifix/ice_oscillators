set term svg
set output 'excitation.svg'
unset key
f(x) = x*x*x/3 - x

set xrange [-2.5 :  2.5]

plot f(x), 'trajectory.dat' using 3:2 with lines lc 0

