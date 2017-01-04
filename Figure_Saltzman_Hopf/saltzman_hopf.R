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
require(glacial)

## glacial package available on request on www.climate.be/itop

t<- seq(-2.0e6,1.5e6,2e3)
tat <- seq(-2.0e6,1.5e6,5e5)/1e3
Rforcing <- sapply(t, function (tt) Insol(ber90(tt),long=90*pi/180,lat=65*pi/180))
Rforcing <- Rforcing - mean(Rforcing)
Rforcing <- approxfun(t,Rforcing)


fmu90 <- function(t) {0.0025 + (t+2.0e6)/2.e6 * (-0.010 )}
## correction 
parSM90$I0=49 
parSM90$sig.mu=0
parSM90$sig.theta=0
parSM90$KR =5.00e-3
## end correction

par <- mapSM90(parSM90)

SM90.Fig1 <- ts_stochastic (ddtSM90,init=c(I=22,mu=253,theta=6.5),t,par=par,
                           Rforcing=zero, stochastic=TRUE, fmu=fmu90)

SM90.Fig2 <- ts_stochastic (ddtSM90,init=c(I=22,mu=253,theta=6.5),t,par=par,
                           Rforcing=Rforcing, stochastic=FALSE, fmu=fmu90)

## correction 
parSM91$I0=49 
parSM91$sig.mu=0
parSM91$sig.theta=0
parSM91$KR =7.530000e-03 
## end correction

fmu91 <- function(t) {0.0015 + (t+1.5e6)/2.e6 * (-0.002 )}
par <- mapSM91(parSM91)

SM91.Fig1 <- ts_stochastic (ddtSM91,init=c(I=22,mu=303,theta=6.5),t,par=par,
                           Rforcing=zero, stochastic=TRUE, fmu=fmu91)


SM91.Fig2 <- ts_stochastic (ddtSM91,init=c(I=22,mu=303,theta=6.5),t,par=par,
                           Rforcing=Rforcing, stochastic=FALSE, fmu=fmu91)


svg('SMtransients.svg', height=4, width=10)
par(mfrow=c(1,2))
plot(t/1.e3, SM90.Fig2[,1], typ='l', frame=FALSE, xlim=c(max(t), min(t))/1e3)
axis(3,at=tat, labels=fmu90(tat*1000), line=2)


plot(t/1.e3, SM91.Fig2[,1], typ='l', frame=FALSE, xlim=c(max(t), min(t))/1e3)
axis(3,at=tat, labels=fmu91(tat*1000), line=2)
dev.off()


