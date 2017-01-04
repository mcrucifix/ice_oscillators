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
BETA  = -1.080 ## asymmetry: Bifurcation for beta = 1 or -1
GAMMA = 0.00 ## forcing amplitude
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

b=c(0.,1.)
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
dev.off ()





