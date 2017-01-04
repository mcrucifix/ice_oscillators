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


PERIOD=0.775701
PP2 = 2*pi/PERIOD
forcing = function(t) sin(PP2*t)
# here the oscillator is programmed to cross a bifurcation. Uncomment second line
# to explore the stationary case. 
Xp = function (t,y,parms) {  
 dx =  - (y[2])   + parms$beta
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

Init = c(0.1,0.1)
times=seq(0,500)*PERIOD/10.
R=seq(length(times)-60, length(times), 1)
## ensemble of experiments varying alpha
parms = list(epsilon = 30.,beta=-0.8)
Sol = lsoda(Init,times,Xp,parms,jacfunc=NULL)
write.table(file="betamo8.txt", Sol[R,], col.names=FALSE, row.names=FALSE)

parms = list(epsilon = 30.,beta=-0.0)
Sol = lsoda(Init,times,Xp,parms,jacfunc=NULL)
write.table(file="betazero.txt", Sol[R,], col.names=FALSE, row.names=FALSE)

parms = list(epsilon = 30.,beta=0.8)
Sol = lsoda(Init,times,Xp,parms,jacfunc=NULL)
write.table(file="betapo8.txt", Sol[R,], col.names=FALSE, row.names=FALSE)


