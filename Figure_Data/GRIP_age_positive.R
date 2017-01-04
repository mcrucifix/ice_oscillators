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

pdf('GRIP-timepositive.pdf', height=6, pointsize=14)

colgrip = seq(1,length(GRIP$V1),10)

GRIP$V1 = - GRIP$V1 
GRIP$V2 =  GRIP$V2 * 1000

plot(GRIP[colgrip,]/1000, typ='l', xlab='time (ka)', ylab='GRIP delta-18 O', frame=FALSE)

dev.off()


