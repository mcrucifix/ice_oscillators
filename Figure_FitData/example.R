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
####################################
# miscellaneous function


plot_trajectory = function(X,plotLR04=TRUE, plotForcing=TRUE, 
                           main='relaxation oscillator as palaeoclimate model',
                           overlay=FALSE)
{
times = X[,1] 
XREF  = X[,2:3]
F     = X[,4]

Xtranslate = function(x) {x*2.}
Ytranslate = function(x) {x-5.0}
Ftranslate = function(x) {x/3.+5.0}
Ltranslate = function(x) {x*2-8}


par(oma=c(0,0,1,1))


plot(times, Xtranslate(XREF[,1]), axes=FALSE, typ='n', lwd=3, ylim=c(-7,6), xlab='Time (ka)', ylab='' )

points(times, Xtranslate(XREF[,1]), type='l')
axis(1, at=(seq(-7,1)*100))

axis(2, at=Xtranslate((-1:1)*1.), lab=c(-1:1)*1., line=1 )


if (plotLR04) 
{
require(glacial)
data(LR04stack)
axis(4, at=Ltranslate(c(3:5)), lab=c(3:5)  ,col='blue')
mtext(side=4, at=1, text = 'LR04 (obs.)', line=2)
points(LR04stack$CE/1000., Ltranslate((LR04stack$delta18O)),col='blue', typ='b', cex=0.6)
}


lines (times, Ytranslate(XREF[,2]), typ='l', lwd=2, col='black')
axis(2, at=Ytranslate(c(-2:2)), labels=c(-2:2), col='black',line=1)

mtext(side=2, at=Xtranslate(0), text='X', line=3, outer=FALSE)
mtext(side=2, at=Ytranslate(0), text='Y', line=3)

if (plotForcing) 
  {
  lines (times, Ftranslate(F), typ='l', lwd=2, col='black')
  labf=c(-2:2)*2
  axis(4, at=Ftranslate(labf), labels=labf, col='black',line=0)
  mtext(side=4, at=Ftranslate(0), text = 'Forcing', line=2)
  }

mtext(main,3,line=3)
}




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
dyn.load('../Fortran/vanderpol.so')
output = as.integer(0)
icalclyap0 = as.integer(0)
icalclyap1 = as.integer(1)



ALPHA = 30.0
BETA  = 0.75
GAMMA = 0.40
RATIO = 2.50

OMEGA_N = 2.20383


PARAMS  = c(ALPHA, BETA, GAMMA,RATIO*OMEGA_N)
TAU = 41 / (2*pi) * (RATIO* OMEGA_N)
na = as.integer(34)

# sets additive fluctuations to zero
b = matrix(c(0.0,0.0), ncol=2, nrow=1) 
ix = as.integer(12)
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

# first : spin-up to find the direction 'ds' of the most unstable perturbation
b=c(0.,0.)
OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)
X = read.table('trajectory.dat')

subsample=seq(1,length(X[,1]),10)
X=X[subsample,]
X[,1] = X[,1] / CONVERT_FACTOR
pdf('example.pdf', height=7, width=6)
plot_trajectory(X)
dev.off()

## ok: now that we have plotted our example, check the number of basins of attraction with this system

###
#source('basin.R')
#B = basin(t0,t1,PARAMS)
## finds two basins

## last point : need to demonstrate the local instability phenomenon 

b=c(0.,0.5)

OUT = .Fortran('propagate_lin', N, state,par,t0,t1, b, ix, output, icalclyap1, ds, lyap,na)

X2 = read.table('trajectory.dat')
X2 = X2[subsample,]
X2[,1] = X2[,1] / CONVERT_FACTOR


pdf('stochastic.pdf', width=6, height=7, pointsize=16)

par(mfrow=c(3,1))
par(mar=c(1,4,1,1))
par(oma=c(4,4,3,0))

plot (X[,1], X[,4], typ='l', lwd=2, col='black', ylab='F(t)',  
                    xlab='', frame=FALSE, axes=FALSE, xlim=c(-700,100))

axis(2)
 

plot(X[,1], X[,2], type='l', frame=FALSE, axes=FALSE, xlim=c(-700,100), 
                   ylim=c(-1.5,1.5),col='black', xlab='', ylab='X', lwd=1.5)

points(X[,1], X2[,2], type='l', frame=FALSE, axes=FALSE, col='red', lwd=1.5)


axis(2, at=(-1:1))
#arrows(-470,-0.4,-470,0, length=0.1, angle=20, col='blue',lwd=1.5)
#arrows(-20,1.4,-20,1.0, length=0.1, angle=20, col='blue',lwd=1.5)


plot(X[,1], X[,3], type='l', frame=FALSE, axes=FALSE, xlim=c(-700,100), 
                         col='black', xlab='', ylab='Y',lwd=1.5)
points(X[,1], X2[,3], type='l', frame=FALSE, axes=FALSE, 
                        xlim=c(-700,0), col='red', lwd=1.5)

axis(1, at=(-7:1)*100)

axis(2)
mtext('time (thousand years)', 1, line=1.5,  outer=TRUE)
dev.off()







