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


#############################################################
# MAIN PROGRAM

# DEFINE MAIN PARAMETERS
# THE LIMIT CYCLE LENGTH OF THE VDP WITH ADOPTED PARAMETERS

# PARAMS : Gamma, Omega
# b : additive noise
# ix: seed
# output : 1 to get the fortran routine to output
#            the trajectory to 'trajectory.dat'; 0 else. 

# natural angular velocity of the oscillator, valid for ALPHA=30.0 and BETA=0.75
dyn.load('vanderpol.so')
output = as.integer(0)
icalclyap0 = as.integer(0)
icalclyap1 = as.integer(1)



ALPHA = 30.0 ## time-scale separation
BETA  = 1.03419 ## asymmetry: Bifurcation for beta = 1 or -1
GAMMA = 0.30 ## forcing amplitude
RATIO = 2.50 ## forcing angular velocity

OMEGA_N = 2.20383 ## constant. Should be calculated for any new value of Beta


PARAMS  = c(ALPHA, BETA, GAMMA,RATIO*OMEGA_N)
TAU = 41 / (2*pi) * (RATIO* OMEGA_N)         ## not sent to the model. Only for diagnostic.
na = as.integer(34)                          ## number of terms for forcing
                                             ## use na=1 for periodic forcing

ix = as.integer(2)
output = as.integer(1)


# define times
#  remember : time units are weird.
#  2*pi = 41,000 years of real time
#  
  CONVERT_FACTOR = 2*pi / 41. 

  t0 = as.numeric(-700.*CONVERT_FACTOR)
  t1 = as.numeric(+100.*CONVERT_FACTOR)

#1. Computes some forward trajector, starting form initial condition 0.2 / -0.4
N = as.integer(1) 
state=matrix(c(0.2, -0.4),N,2)

par = t(matrix(rep(PARAMS,N), ncol=4, nrow=1))

output = as.integer(1)

# make sure that par, lyap and ds have the right form.
lyap  =  state[,1]
ds    =  state
lyap[] = 0
ds  [,1] = 1.
ds  [,2] = 0.

# compute four sample trajectories with different parameters

b=c(0.,0.)
par[2] = 0.75
par[3] = 0.40
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
X_B0_BETA075 = read.table('trajectory.dat')
subsample=seq(1,length(X_B0_BETA075[,1]),10)
X_B0_BETA075 = X_B0_BETA075[subsample,]

b=c(0.,1.00)/sqrt(RATIO*OMEGA_N)
ix=12
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
X_B1_BETA075 = read.table('trajectory.dat')[subsample,]

b=c(0.,0.)
par[2] = 1.035
par[3] = 0.
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
X_B0_BETA103 = read.table('trajectory.dat')[subsample,]

b=c(0.,0.15)/sqrt(par[4])
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
X_B1_BETA103 = read.table('trajectory.dat')[subsample,]


b=c(0.,0.30)/sqrt(par[4])
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
X_B2_BETA103 = read.table('trajectory.dat')[subsample,]

## a few additional trials to check the formula published in the paper (hand edits)
par[1] = 30
par[2] = 1.07
b=c(0.,0.010)/sqrt(par[4])
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
X_B3_BETA103 = read.table('trajectory.dat')[subsample,]


### note : 
## ATTENTION !!! 
## about the actual value of b. 

## in the code, the equation is code as follows:
## dy = f(x,y)/OMEGA  dt + F(t) + b * dw, then
## in /// adimensional // time units.
## While in the paper, the equation is as follows:
## dy = f(x,y)  dt + F(OMEGA * t) + b' * dw, then
## we get : b' = b * SQRT(OMEGA).
## This, the 'b'' of the paper should be first divided
## by SQRT(OMEGA) before being passed on to the code

## If now we write the DIMENSIONAL equations
## say as :

## dy = 1/TAU * [ f(x,y)  dt' + F(OMEGA * t' / TAU )] + 1/SQRT(TAU) b' * dw', then

## then the equation is equivalent as above, i. e., setting t' = t * TAU, 
## i. e. : dt' = dt * TAU and dw' = dw * sqrt(TAU) 
## we have
## dy =  [ f(x,y)  dt' + F(OMEGA * t' )] +  b' * dw', then


postscript('local_instability.ps', height=6, width=5)

TIMES = X_B0_BETA103[,1] / CONVERT_FACTOR

plot(TIMES, X_B0_BETA075[,2], typ='l', lty=1, lwd=2, 
           col='black', axes=F, ylim=c(-5,5), 
           xlab='time (thousand years)', ylab='',
           xlim=c(-700,100), yaxs='i')
lines(TIMES, X_B0_BETA075[,3]/2.-3, lty=1, lwd=2, col='black')
lines(TIMES, X_B0_BETA075[,4]*GAMMA+3, lty=1, lwd=2, col='black')

lines(TIMES, X_B1_BETA075[,2], lty=2, lwd=2,   col='red')
lines(TIMES, X_B1_BETA075[,3]/2.-3, lty=2, lwd=2, col='red')

axis(1, at=seq(-7,1)*100)

axis(2, at=seq(-1,1))
axis(2, at=seq(-1,1)-3., labels=seq(-1,1)/2.)
axis(2, at=seq(-1,1)+3., labels=seq(-1.,1.))

mtext('x', 2, at=0, line=3)
mtext('y', 2, at=-3, line=3)
mtext('F(t)', 2, at=+3, line=3)

dev.off() 


postscript('excitability.ps', height=6, width=5)

TIMES = X_B0_BETA103[,1] 

plot(TIMES, X_B0_BETA103[,2], typ='l', lty=1, lwd=2, 
           col='black', axes=F, ylim=c(-5,2), 
           xlab='time (arbitrary_units)', ylab='',
           yaxs='i')
lines(TIMES, X_B0_BETA103[,3]/2.-3, lty=1, lwd=2, col='black')

lines(TIMES, X_B1_BETA103[,2],      lty=2, lwd=2,   col='red')
lines(TIMES, X_B1_BETA103[,3]/2.-3, lty=2, lwd=2,   col='red')

lines(TIMES, X_B2_BETA103[,2],      lty=3, lwd=2,   col='blue')
lines(TIMES, X_B2_BETA103[,3]/2.-3, lty=3, lwd=2,   col='blue')


axis(1)

axis(2, at=seq(-1,1))
axis(2, at=seq(-1,1)-3., labels=seq(-1,1)/2.)

mtext('x', 2, at=0, line=3)
mtext('y', 2, at=-3, line=3)


legend('topright', c('b=0','b=0.15','b=0.30'), col=c('black','red','blue'), lty=1)
# mtext('Effect of stochastic perturbations on excitable system', 3, line=2)
#mtext('b=0.35', font=13, 3, line=1)

dev.off ()





