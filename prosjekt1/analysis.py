import numpy as np
import matplotlib.pyplot as plt

# Read computed data set
fil = "n10000.txt"
data_file = open(fil, 'r')
lines = data_file.readlines()

n = len(lines)
v_0 = v_n_1 = 0

v = [v_0]
for line in lines:
    v.append(float(line))

v.append(v_n_1)
v = np.array(v)

# Compute analytical solution
x_0 = 0.
x_n_1 = 1.
h = (x_n_1-x_0)/(n+1)
x = []
# Compute x-values in same manner as in the c++ program
# so as to keep a potential error here constant
for i in range(n+2):
    x.append(i*h)
x = np.array(x)

def u_analytical(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

u = u_analytical(x)

eps = 0
for i in range(len(u)):
    if u[i]> 10**(-15):
        new_eps = np.log10(abs((v[i] - u[i])/u[i]))
        if abs(new_eps)>eps:
            eps = new_eps
print eps

"""
plt.plot(x, v, x, u)
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend("Numerical, Analytical")
plt.title("n=1000")
plt.show()
"""
