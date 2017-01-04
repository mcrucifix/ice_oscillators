## ---------------------------------------------------------
## Code realeased under Creative Common Licence
## The code was used to generate the figures available in the Paper:
## "Oscillators and relaxation phenomena in Pleistocene climate theory "
## published in the Philosophical Transactions of the Royal Society
## Code delivered without Garantee and Support whatsoever
## Contact michel.crucifix for details 
## Published on arXiv.org on 21 DEC 2012.
## Michel Crucifix
## I keep intellectual ownership of the work
## ---------------------------------------------------------
## modified van der pol oscillator to study glacial cycles

phi = function (x) {x^3/3-x}
phi_prime = function (x) {x^2-1}


## van der pol with three stable states and two unstable states 
require(odesolve)
require(glacial)

 # defines phi function
# phi = function (x) { 
# (x+1.7)*(x+1.58)*(x+0.8)*(x-0)*(x-0.50)}


#forcing = function(t) {o = ber90(t*10000) 
#                       100*o["ecc"]*sin(o["varpi"]) + 100 * (o["eps"]-0.409)} 


PERIOD=0.9
PP2 = 2*pi/PERIOD
forcing = function(t) sin(PP2*t)
# here the oscillator is programmed to cross a bifurcation. Uncomment second line
# to explore the stationary case. 
Xp = function (t,y,parms) {  
 dx =  - (y[2])   + 0.50* forcing(t)
 dy =  - parms$epsilon*(phi(y[2]) - y[1])   
# dx =  -parms$epsilon * (y[2]-0.00+parms$alpha+(t+50)/200)   + 0.15 * forcing(t)
 list(c(dx,dy))
}

# Jacobian, to speed up integration in lsoda. 
Jac = function (t,y,parms) {
 dxdx =  0
 dxdy =  -1
 dydy =  - parms$epsilon * ( phi_prime(y[2]))
 dydx =  1
 c(dxdx,dxdy,dydx,dydy)
}


LENGTH=20

Init = c(0.1,0.1)
times=seq(0,800)*PERIOD/LENGTH

## ensemble of experiments varying alpha
parms = list(epsilon = 11.1, alpha=0)
Sol = lsoda(Init,times,Xp,parms,jacfunc=NULL)


Seq1= 67 + seq(1,5)*LENGTH*3
Seq2= Seq1 + LENGTH
Seq3= Seq2 + LENGTH
Seq4= c(Seq1, Seq2, Seq3)

pdf('Locking.pdf', 5, 6, pointsize=11)
plot(times, Sol[,2], type='l', frame=FALSE, axes=FALSE, ylim=c(-8,5), xlim=c(5,12), lwd=2, xlab='time (Arb. Units)', ylab='')
lines(times, forcing(times)+3, type='l')
lines(times+PERIOD,   Sol[,2]-3, type='l', lwd=2)
lines(times+2*PERIOD, Sol[,2]-6, type='l', lwd=2)

axis(1)


for (i in Seq1) {lines(times[c(i,i)], c(Sol[i,2], forcing(times[i])+3), lty=2)
                 points(times[i], forcing(times[i])+3, pch=16, col='yellow')
                 points(times[i], forcing(times[i])+3, pch=1)
                 points(times[i], Sol[i,2], pch=16, col='yellow')
                 points(times[i], Sol[i,2], pch=1)}
for (i in Seq2) {lines(times[c(i,i)], c(Sol[i-LENGTH,2]-3, forcing(times[i])+3), lty=2)
                 points(times[i], forcing(times[i])+3, pch=17, col='blue')
                 points(times[i], forcing(times[i])+3, pch=2)
                 points(times[i], Sol[i-LENGTH,2]-3, pch=17, col='blue')
                 points(times[i], Sol[i-LENGTH,2]-3, pch=2)}
for (i in Seq3) {lines(times[c(i,i)], c(Sol[i-LENGTH-LENGTH,2]-6, forcing(times[i])+3), lty=2)
                 points(times[i], forcing(times[i])+3, pch=18, col='red')
                 points(times[i], forcing(times[i])+3, pch=5)
                 points(times[i], Sol[i-LENGTH-LENGTH,2]-6, pch=18, col='red')
                 points(times[i], Sol[i-LENGTH-LENGTH,2]-6, pch=5)}

text(5,4.5, 'Forcing', pos=4)
text(5,1.5, 'x(t) : Locking 1', pos=4)
text(5,-1.5, 'x(t) : Locking 2', pos=4)
text(5,-4.5, 'x(t) : Locking 3', pos=4)

dev.off()
