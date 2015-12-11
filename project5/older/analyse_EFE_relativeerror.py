from numpy import *
from matplotlib.pyplot import *

infile = open("EFET0p5Nt5000Nx200.txt", 'r')

parameters = infile.readline()
v_num = []

for line in infile:
    v_num.append(float(line))

x = linspace(0, 1, 201)

def v_an(x):
    T=0.5
    v = 0
    for n in range(1,200):
        v += (1./n)*sin(n*pi*x)*exp(-T*(n*pi)**2)
    v *= -2/pi
    return v


u_num=array(v_num) + 1 - x
u_an=array(v_an(x)) + 1 - x
print u_an[0:4]
reler = abs((u_an-u_num))
#print reler[0]
print len(reler)

"""
#plot(x,u_num, x, u_an)
plot(x[:-1], 1E6*abs((u_an-u_num)[:-1]/u_an[:-1]))
#legend(["Numerical", "Analytic"])

override = {
    'fontsize'            : 'large',
    'verticalalignment'   : 'baseline',
    'horizontalalignment' : 'center'
    }
xlabel("x", override)
ylabel("Relative error $\\times$1E-6", override)
title("Explicit Forward Euler Relative Error \n T=0.5 Nt=5000 Nx=200")
show()
"""
"""
ylabel("v(x, t=5")
title("Explicit Forward Euler \n Nx=20 Nt=5000")
show()
"""
