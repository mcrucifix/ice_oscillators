set term postscript enhanced color eps
set output "saltzman.eps"

coef=0.2
v1=0.8
s=0.8
r=0.6
p=1.
phi(y)=y*y*y+s*y*y-r*y
dx1(x,y)=-x*v1 - y
dy1(x,y)=-(phi(y)-p*x)



set parametric
set table "p1.tmp"; plot [-2:2] phi(t),t; unset table
set table "p2.tmp"; plot [-2:2] t,-t*v1; unset table
set isosam 11,11
unset parametric
set xrange [-1.5:1.5]
set yrange [-2:2]
set table "field2xy.tmp"
splot x,y w l
unset table

set samples 20
#set label '$y=\phi(x)$' at 0.38,1.69 textcolor rgb "brown"
#set label '$x*(1-v)=y$' at -0.76,0.38 textcolor rgb "blue"
#set label '$x$' at 1.8,0.1
#set label '$y$' at 0.1,2.
set label 'XPHI' at 0.38,1.69 textcolor rgb "brown"
set label 'X' at 1.4,0.1
set label 'Y' at 0.1,2.
set arrow from 0,-2 to 0,2 lt rgb "black"
set arrow from -1.5,0 to 1.5,0 lt rgb "black"
unset xzeroaxis
unset yzeroaxis
unset border
unset xtics
unset ytics
unset border
unset key
set parametric
set trange [-2:2]
plot  "field2xy.tmp" u 1:2:(1./7*coef*dx1($1,$2)):(1./7*coef*dy1($1,$2)) w vec lt rgb "grey",  \
       "p1.tmp" with lines lw 4 lt 1 lc rgb "brown", \
       "p2.tmp" with lines lw 4 lt 1 lc rgb "blue"
