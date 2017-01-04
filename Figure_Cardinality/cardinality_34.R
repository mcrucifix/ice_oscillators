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

# miscellaneous function

#dyn.unload('vanderpol.so')
dyn.load('../Fortran/vanderpol.so')
source('../src/cluster.R')
source('../src/basin.R')

output = as.integer(0)
icalclyap0 = as.integer(0)
icalclyap1 = as.integer(1)

OMEGA_1 = 2.20539 
ALPHA = 11.111
BETA = 0.25
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


outputdata=sprintf('Cardinality_%i.RData',na)
outputpng =sprintf('Cardinality_%i.png',na)

if (!file.exists(outputdata)) {local({

# print('Data file not in current directory; need to recreate it')

tI = as.numeric(-40.*2.*pi)
t0 = as.numeric(0.0)

Gamma=seq(0.1,4.1,0.01)
Omega=seq(0.5,3.5,0.01)*OMEGA_1

N=length(Gamma)
M=length(Omega)

K=matrix(NA, ncol=N, nrow=M)
 for (i in (1:N)) for (j in (1:M)) {print(i)
       K[j,i]=basin(tI,t0, c(ALPHA, BETA, Gamma[i], Omega[j]))$nc}

KList=list(x=Omega/OMEGA_1, y=Gamma, z=K)

save(KList, file=outputdata)
})}


load(outputdata)
png (outputpng,width=5,height=5,res=600,unit='in')

Col=c('grey','blue','red','green', 'gold','white')
breaks=(0:6)
image(KList, breaks=breaks, col=Col, xlab='', ylab=expression(gamma), 
main='Cardinality of pullback attractor : astro forcing', frame.plot=TRUE, axes=FALSE, xlim=c(0.5,3.5))
axis(2)
axis(1,at=(1:4),lab=(1:4))
legend('topleft', legend=c('1','2','3','4','5','>5'), bg='white', col=c(Col[-6],'black'), fill=Col, text.col='black')

mgps = c(3,0.2,0)

p1 = 23716/41000
p2 = 22428/41000
p3 = 18976/41000

mtext(expression(omega[epsilon*1]/omega[N]), line=2, col='black',1)

axis(1, at=(1:5)*p1, lab=(1:5), col='blue', col.tick='blue',col.axis='blue',cex.lab=0.6, cex.axis=0.6, col.lab='blue', tcl=-0.2,mgp=mgps)
mtext(expression(omega[p*1]/omega[N]), line=3, col='blue',1, at=1.2)


axis(1, at=(2:5)*p2, lab=(2:5), col='red', col.tick='red',col.axis='red',cex.lab=0.6, cex.axis=0.6, col.lab='blue', tcl=-0.2,mgp=mgps)
mtext(expression(omega[p*2]/omega[N]), line=3, col='red', 1, at=2.0)


axis(1, at=(2:4)*p3, lab=(2:4), col='purple', col.tick='purple',col.axis='purple',cex.lab=0.6, cex.axis=0.6, col.lab='blue', tcl=-0.2,mgp=mgps)
mtext(expression(omega[p*3]/omega[N]), line=3, col='purple', 1, at=2.8)



dev.off()

