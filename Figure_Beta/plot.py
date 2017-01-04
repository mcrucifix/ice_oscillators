#!/usr/bin/env python

from pylab import *
from numpy import *

params = { 'figure.figsize': [5,4]}
rcParams.update(params)



x,y = loadtxt("Beta_Period.txt", unpack=True)
t,dummy,s1 = loadtxt("betamo8.txt", unpack=True)
t,dummy,s2 = loadtxt("betazero.txt", unpack=True)
t,dummy,s3 = loadtxt("betapo8.txt", unpack=True)

# the main axes is subplot(111) by default
plot(x, y, lw=2)
axis([-1, 1,0,7])
xlabel('$\\beta$')
ylabel('Natural Oscillating Period (arbitrary units)')
title('')

# this is an inset axes over the main axes
a = axes([.2, .6, .2, .2], axisbg=None)
plot(t, s1)
title('$\\beta = -0.8$')
xlabel('$t$')
ylabel('$y$')
setp(a, xticks=[], yticks=[])

# this is another inset axes over the main axes
a = axes([0.6, 0.6, .2, .2], axisbg=None)
plot(t, s3)
title('$\\beta = 0.8$')
xlabel('$t$')
ylabel('$y$')
setp(a, xticks=[], yticks=[])

a = axes([0.4, 0.16, .2, .2], axisbg=None)
plot(t, s2)
title('$\\beta = 0.0$')
xlabel('$t$')
ylabel('$y$')
setp(a, xticks=[], yticks=[])

savefig("beta.pdf", transparent=True)
