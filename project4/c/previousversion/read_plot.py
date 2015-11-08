import numpy as np
import matplotlib.pyplot as plt

f = "c_ordered_L20temp24mcE2toE7.txt"
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


# Mean energy
plt.figure(1)
plt.subplot(211)
plt.plot(np.log10(number_mcs), E_exp, 'o-')
plt.ylabel("Energy/J")
plt.title("Expectation values approaching equilibrium\n L=20, T=2.4, init. config. = ordered")
plt.legend(["<E>"])
plt.grid()


# Mean magnetization
plt.subplot(212)
plt.plot(np.log10(number_mcs), M_exp, 'o-')
plt.plot(np.log10(number_mcs), M_abs_exp, 'o-')
plt.ylabel("Magnetization")
plt.legend(["|<M>|","<|M|>"], loc=4)
plt.grid()

plt.xlabel("log10(Time [number of MC cycles])")
plt.show()
