import numpy as np
import matplotlib.pyplot as plt

f = "test2.txt"
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
accepted_counter = []

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
    chi_abs.append(float(line[9]))
    accepted_counter.append(float(line[10]))  # includes \n. Python handles it


n_mcs = np.array(number_mcs[0:6])
count1 = np.array(accepted_counter[0:6])/n_mcs
count2 = np.array(accepted_counter[6:12])/n_mcs
count3 = np.array(accepted_counter[12:])/n_mcs
n_mcs = np.log10(n_mcs)


plt.figure(1)
plt.subplot(311)
plt.plot(n_mcs, count3, 'ro-')
plt.title("Number of accepted spin flips\nL=20")
plt.legend(["T=3\ninit. config.=random"], loc=4)
plt.grid()

plt.subplot(312)
plt.plot(n_mcs, count2, 'go-')
plt.ylabel("Number of accepted moves")
plt.legend(["T=2\ninit. config.=random"], loc=4)
plt.grid()

plt.subplot(313)
plt.plot(n_mcs, count1, 'bo-')
plt.legend(["T=1\ninit. config.=ordered"], loc=1)
plt.grid()

plt.xlabel("log10(Time [number of MC cycles])")
plt.show()
