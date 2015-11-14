import numpy as np
import matplotlib.pyplot as plt

f = "e_L80dt03mcE6proc0.txt"
data_file = open(f, 'r')

L          = []
number_mcs = []
time       = []
temp       = []
E_exp      = []
M_exp      = [] # |<M>|
M_abs_exp  = [] # <|M|>
Cv         = []
chi        = []
chi_abs    = []

for line in data_file:
    line = line.split()
    L.append(float(line[0]))
    number_mcs.append(float(line[1]))
    time.append(float(line[2]))
    temp.append(float(line[3]))
    E_exp.append(float(line[4]))
    M_exp.append(float(line[5]))
    M_abs_exp.append(float(line[6]))
    Cv.append(float(line[7]))
    chi.append(float(line[8]))
    chi_abs.append(float(line[9])) # includes \n. Python handles it

f = "e_L80dt03mcE6proc1.txt"
data_file = open(f, 'r')

for line in data_file:
    line = line.split()
    L.append(float(line[0]))
    number_mcs.append(float(line[1]))
    time.append(float(line[2]))
    temp.append(float(line[3]))
    E_exp.append(float(line[4]))
    M_exp.append(float(line[5]))
    M_abs_exp.append(float(line[6]))
    Cv.append(float(line[7]))
    chi.append(float(line[8]))
    chi_abs.append(float(line[9])) # includes \n. Python handles it


plt.plot(temp, E_exp)
plt.xlabel("$\\alpha$")
plt.show()

plt.plot(temp, M_abs_exp)
plt.show()

plt.figure(1)
plt.subplot(211)
plt.plot(temp, Cv)
plt.xlabel("T")
plt.ylabel("$C_v$")
plt.title("L=60, number of cycles = 1E5")

plt.subplot(212)
plt.plot(temp, chi_abs)
plt.xlabel("T")
plt.ylabel("$\\chi$")
plt.show()


"""
# Mean energy
plt.figure(1)
plt.subplot(211)
plt.plot(number_mcs, E_exp, '-')
plt.ylabel("Energy/J")
plt.title("Expectation values approaching equilibrium\n L=20, T=2.4, init. config. = random")
plt.legend(["<E>"], loc=3)
plt.grid()

# Mean magnetization
plt.subplot(212)
plt.plot(number_mcs, M_exp, '-')
plt.plot(number_mcs, M_abs_exp, '-')
plt.ylabel("Magnetization")
plt.legend(["|<M>|","<|M|>"], loc=2)
plt.grid()

plt.xlabel("log10(Time [number of MC cycles])")
plt.show()
"""


