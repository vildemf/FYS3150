from matplotlib.pyplot import *
from numpy import *

alpha = 2.
r = linspace(0, 4, 10000)
psi = exp(-alpha*r)

plot(r, psi)
xlabel("$r$")
ylabel("$\psi_{1s}(r)$")
title("Single particle wave function")
show()
