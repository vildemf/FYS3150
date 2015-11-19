from numpy import *
from matplotlib.pyplot import *

def v(x):
    T=5
    v = 0
    for n in range(1,10):
        v += (1./n)*sin(n*pi*x)*exp(-T*(n*pi)**2)

    v *= -2*pi
    return v

x = linspace(0,1,21)



plot(x, v(x))
show()
