set term postscript enhanced color eps
set output "vdp.eps"
phi(y) = (y)**3/3 - (y)
set parametric
set table "p1.tmp"; plot [-2:-1] phi(t),t; unset table
set table "p2.tmp"; plot [-1:1] phi(t),t; unset table
set table "p3.tmp"; plot [1:2] phi(t),t; unset table
set isosam 11,11
unset parametric
set xrange [-1.5:1.5]
set yrange [-2:2]
set table "field2xy.tmp"
splot x,y w l
unset table

coef=0.2
beta=-0.80
F=0.6
dx1(x,y)=-(y+beta)
dy1(x,y)=-(phi(y)-x)

set samples 20
#set label '$y=\phi(x)$' at 0.38,1.69 textcolor rgb "brown"
#set label '$x=\beta$' at -0.76,0.38 textcolor rgb "blue"
#set label '$x$' at 1.8,0.1
#set label '$y$' at 0.1,2.
set label 'XPHI' at 0.42,1.69 textcolor rgb "brown"
set label 'YBETA' at -0.56,-beta+0.1 textcolor rgb "blue"
set label 'YBETAFM' at 0.46,-beta+F/3 textcolor rgb "blue"
set label 'X' at 1.4,0.1
set label 'Y' at 0.1,2.
set arrow from 0,-2 to 0,2 lt rgb "black"
set arrow from -1.5,0 to 1.5,0 lt rgb "black"
set arrow from 0.4,-beta to 0.4,-beta+F lt  rgb "blue"
set arrow from 0.4,-beta to 0.4,-beta-F lt  rgb "blue"
unset xzeroaxis
unset yzeroaxis
unset border
unset xtics
set yzeroaxis
set ytics -1, 1, 1
unset border
unset key
set parametric
set trange [-2:2]
plot  "field2xy.tmp" u 1:2:(1./10*coef*dx1($1,$2)):(coef*dy1($1,$2)) w vec lt rgb "grey",  \
       "p1.tmp" with lines lw 4 lt 1 lc rgb "brown", \
       "p2.tmp" with lines lw 4 lt 2 lc rgb "brown", \
       "p3.tmp" with lines lw 4 lt 1 lc rgb "brown", \
       t,-beta lw 4  lt 1 lc rgb "blue", \
       t,-beta+F lw 4  lt 2 lc rgb "blue", \
       t,-beta-F lw 4  lt 2 lc rgb "blue", \
       phi(-beta),-beta with points pt 6 ps 5 lc rgb "red"
