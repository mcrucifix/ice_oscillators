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


#forcing
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


PERIOD=0.775701
PP2 = 2*pi/PERIOD
forcing = function(t) sin(PP2*t)
# here the oscillator is programmed to cross a bifurcation. Uncomment second line
# to explore the stationary case. 
Xp = function (t,y,parms) {  
 dx =  - ( y[2])  # + 0.50* forcing(t)
 dy =  - parms$epsilon*(phi(y[2]) -   y[1] )  
 dz =  -  5*(y[3] + y[1] - 0.0*y[2] )
# dx =  -parms$epsilon * (y[2]-0.00+parms$alpha+(t+50)/200)   + 0.15 * forcing(t)
 list(c(dx,dy,dz))
}

# Jacobian, to speed up integration in lsoda. 
Jac = function (t,y,parms) {
 dxdx =  0.
 dxdz =  0.
 dxdx =  1
 dydx =  1
 dydy =  - parms$epsilon * ( phi_prime(y[2]))
 dydz =  0.
 dzdx =  1.
 dzdy =  0.
 dzdz =  -1.

 c(dxdx,dxdy,dxdz,dydx,dydy, dydz, dzdx, dzdy, dzdz)
}

Init = c(0.1,0.1, 0.1)
times=seq(0,800)*PERIOD/40.

## ensemble of experiments varying alpha
parms = list(epsilon = 30., alpha=0)
Sol = lsoda(Init,times,Xp,parms,jacfunc=NULL)


PERIOD=0.775701
PP2 = 2*pi/PERIOD
Xp = function (t,y,parms) {  
 dx =  - ( y[2] ) #  + 0.50* forcing(t)
 dy =  - parms$epsilon*(phi(0.5 * y[2] + 0.5 * y[3]) - y[1])  
 dz =  -  10*(y[3] - y[2]  )

# dx =  -parms$epsilon * (y[2]-0.00+parms$alpha+(t+50)/200)   + 0.15 * forcing(t)
 list(c(dx,dy,dz))
}


Init = c(0.1,0.1, 0.1)
times=seq(0,800)*PERIOD/40.

## ensemble of experiments varying alpha
parms = list(epsilon = 30., alpha=0)
Sol = lsoda(Init,times,Xp,parms,jacfunc=NULL)

pdf('Overshoot.pdf', 5, 6, pointsize=11)
plot(times, Sol[,2], type='l', frame=FALSE, axes=FALSE, ylim=c(-14,3), xlim=c(5,15), lwd=2, xlab='times (Arb. Units)', ylab='', title='nominal model with overshoot dynamics')
lines(times, Sol[,2], type='l')
lines(times, Sol[,3]-6, type='l', lwd=2)
lines(times, Sol[,4]-12, type='l', lwd=2)

mtext('x', side=2, at=0, line=3)
mtext('y', side=2, at=-6, line=3)
mtext('z', side=2, at=-12, line=3)

axis(2, at=seq(-2,2),    labels=seq(-2,2))
axis(2, at=seq(-2,2)-6,  labels=seq(-2,2))
axis(2, at=seq(-2,2)-12, labels=seq(-2,2))

axis(1)



dev.off()
