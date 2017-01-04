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

#dyn.unload('vanderpol.so')

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


outputdata=sprintf('Cardinality_beta06_alpha_11.RData',na)
outputpng =sprintf('Cardinality_beta96_alpha_11.png',na)

if (!file.exists(outputdata)) {local({

dyn.load('../Fortran/vanderpol.so')
source('../src/cluster.R')
source('../src/basin.R')
# print('Data file not in current directory; need to recreate it')

tI = as.numeric(-40.*2.*pi)
t0 = as.numeric(0.0)

Gamma=seq(0.1,1.0,0.005)
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

par(mar=c(8,4,3,2))

Col=c('grey','blue','red','green', 'gold','white')
breaks=(0:6)
image(KList, breaks=breaks, col=Col, xlab='', ylab=expression(gamma), 
main='Number of possible solutions', frame.plot=TRUE, axes=FALSE)
mtext('astronomical forcing  - beta=0.60', 3, line=0)
points(c(2.525), c(0.35) , pch=4, cex=1, col='black')
axis(2)
legend('topleft', legend=c('1','2','3','4','5','>5'), bg='white', col=c(Col[-6],'black'), fill=Col, text.col='black')



mgps = c(3,0.2,0)

p1 = 23716/41000
p2 = 22428/41000
p3 = 18976/41000



axis(1,at=((0:6)*25/ 41),lab=(0:6)*25,line=0, padj=-0.5)     # main axis, in real time

mtext(expression(ka), line=0.6,1, at=3.3, cex=1)

mtext('natural period', 1, line=6)

axis(1, at=(1:3), lab=(1:3), col='red', col.tick='red',
        col.axis='red',cex.lab=0.6, cex.axis=0.6, 
        col.lab='red', tcl=-0.2, mgp=mgps,line=2, padj=-0.65) 
mtext(expression(p[epsilon*1]), line=2, col='red',1, at=3.3, cex=0.6)

axis(1, at=(1:5)*p1, lab=(1:5), col='blue', col.tick='blue',
        col.axis='blue',cex.lab=0.6, cex.axis=0.6, 
        col.lab='blue', tcl=-0.2, mgp=mgps,line=3,padj=-0.65) 

mtext(expression(p[p*1]), line=3, col='blue',1, cex=0.6, at=5.2*p1)

axis(1, at=(1:5)*p2, lab=(1:5), col='blue', col.tick='blue',
        col.axis='blue',cex.lab=0.6, cex.axis=0.6, 
        col.lab='blue', tcl=-0.2, mgp=mgps,line=4, padj=-0.65) 

mtext(expression(p[p*2]), line=4, col='blue',1, cex=0.6, at=5.2*p2)

axis(1, at=(1:5)*p3, lab=(1:5), col='blue', col.tick='blue',
        col.axis='blue',cex.lab=0.6, cex.axis=0.6, 
        col.lab='blue', tcl=-0.2, mgp=mgps,line=5,padj=-0.65) 

mtext(expression(p[p*3]), line=5, col='blue',1, cex=0.6, at=5.2*p3)



dev.off()

