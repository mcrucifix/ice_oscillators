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
GRIP = read.table('ngrip-d18o-50yr.txt', skip=80)
GRIP.colnames=c('age (yr)', 'd18O')

require(glacial)
data(LR04stack)

svg('GRIP-LR04.svg')

par(mfrow=c(2,1))

colgrip = seq(1,length(GRIP$V1),10)

plot(GRIP[colgrip,]/1000, typ='l', xlab='age (ka)', ylab='GRIP delta-18 O', frame=FALSE)

plot(-LR04stack$CE/1000., LR04stack$delta18O, xlim=c(0., 1000), typ='l', xlab='age (ka)', ylab='LR04 stack', frame=FALSE)

dev.off()


