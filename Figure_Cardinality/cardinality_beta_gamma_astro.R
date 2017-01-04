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


output = as.integer(0)
icalclyap0 = as.integer(0)
icalclyap1 = as.integer(1)

OMEGA_1 = 2*pi/3.18586
ALPHA = 11.1
BETA = 0.60
GAMMA = 0.5
# sets additive fluctuations to zero
b = matrix(c(0.0,0.0), ncol=2, nrow=1) 
ix = as.integer(1)

# define times
#  remember : time units are weird.
#  2*pi = 41,000 years of real time
#  
 CONVERT_FACTOR = 2*pi / 41. 
# na : 1 for periodic forcing ; 34 for astronomical forcing
na=as.integer(34)


outputdata=sprintf('Cardinality_beta_gamma_astro.RData',na)
outputpng =sprintf('Cardinality_beta_gamma_astro.png',na)

if (!file.exists(outputdata)) {local({

#dyn.unload('vanderpol.so')
dyn.load('../Fortran/vanderpol.so')
source('../src/cluster.R')
source('../src/basin.R')
# print('Data file not in current directory; need to recreate it')

tI = as.numeric(-40.*2.*pi)
t0 = as.numeric(0.0)

Gamma=seq(0.01,1.01,0.005)
Beta=seq(-1.5,1.5,0.01)
Omega=2.5*OMEGA_1

N=length(Gamma)
M=length(Beta)

K=matrix(NA, ncol=N, nrow=M)
 for (i in (1:N)) for (j in (1:M)) {print(i)
       K[j,i]=basin(tI,t0, c(ALPHA, Beta[j], Gamma[i], Omega))$nc}

KList=list(x=Beta, y=Gamma, z=K)

save(KList, file=outputdata)
})}


load(outputdata)
png (outputpng,width=5,height=5,res=600,unit='in')

par(mar=c(4,4,3,2))

Col=c('grey','blue','red','green', 'gold','white')
breaks=(0:6)
image(KList, breaks=breaks, col=Col, xlab=expression(beta), ylab=expression(gamma), 
main='astronomical forcing', frame.plot=TRUE, axes=FALSE)
#mtext('astronomical forcing  ', 3, line=0)
axis(2)
legend('topleft', legend=c('1','2','3','4','5','>5'), bg='white', col=c(Col[-6],'black'), fill=Col, text.col='black')



mgps = c(3,0.2,0)

axis(1)

mtext(expression(ka), line=0.6,1, at=3.3, cex=1)

mtext(expression(beta), 1, line=6)



dev.off()

