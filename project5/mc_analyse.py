from numpy import *
from matplotlib.pyplot import *

filename = open("testing.txt", "r")

average = []
for line in filename:
    #line = line.split()
    average.append(float(line))

x = linspace(0, 1, 226)
plot(x, average)
#axis([0, 1, ])
show()