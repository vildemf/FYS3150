import numpy as np
import matplotlib.pyplot as plt

# Read computed data set
def read_psi(filename):
    fil = filename
    data_file = open(fil, 'r')
    lines = data_file.readlines()
    n = len(lines)
    psi = [0] # psi(0) = 0
    for line in lines:
        psi.append(float(line))
    psi = np.array(psi)
    return psi, n

psi0, n = read_psi("E0.txt")
psi1 = read_psi("E1.txt")[0]
psi2 = read_psi("E2.txt")[0]

rho_min = 0
rho_max = 8
N_step = n+2
h = (rho_max - rho_min)/float(N_step)
rho = np.linspace(rho_min, rho_max - h, n+1)

plt.plot(rho, abs(psi0)**2, rho, abs(psi1)**2, rho, abs(psi2)**2)
plt.xlabel("rho")
plt.ylabel("|psi|**2")
plt.legend(["n=0", "n=1", "n=2"])
plt.title("Exercise d) omega = 1.5, w/Coulomb  (N_step=200, rho_max=8)")
plt.show()
