import numpy as np
import matplotlib.pyplot as plt

def read_2cpus(f):
    f1=f[0]
    f2=f[1]
    data_file = open(f1, 'r')
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
        
    data_file = open(f2, 'r')
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
    return temp, E_exp, M_abs_exp, Cv, chi_abs

filenames = [("e_L20dt02mcE6proc0.txt", "e_L20dt02mcE6proc1.txt"), ("e_L40dt02mcE6proc0.txt", "e_L40dt02mcE6proc1.txt"),\
             ("e_L60dt03mcE6proc0.txt", "e_L60dt03mcE6proc1.txt"), ("e_L80dt03mcE6proc0.txt", "e_L80dt03mcE6proc1.txt")]

temp20, E20, M20, C20, X20 = read_2cpus(filenames[0])
temp40, E40, M40, C40, X40 = read_2cpus(filenames[1])
temp60, E60, M60, C60, X60 = read_2cpus(filenames[2])
temp80, E80, M80, C80, X80 = read_2cpus(filenames[3])

plt.figure(1)
plt.subplot(211)
plt.plot(temp20, E20)
plt.plot(temp40, E40)
plt.plot(temp60, E60)
plt.plot(temp80, E80)
plt.ylabel("$\\langle E\\rangle$/J")
plt.title("Model of Phase Transition \n Number of MC cycles = 1E6")
plt.legend(["L=20", "L=40", "L=60", "L=80"], loc=2)

plt.subplot(212)
plt.plot(temp20, M20)
plt.plot(temp40, M40)
plt.plot(temp60, M60)
plt.plot(temp80, M80)
plt.xlabel("Temperature Tk/J")
plt.ylabel("$\\langle |M| \\rangle$")
plt.legend(["L=20", "L=40", "L=60", "L=80"], loc=3)
plt.show()

plt.figure(2)
plt.subplot(211)
plt.plot(temp20, C20)
plt.plot(temp40, C40)
plt.plot(temp60, C60)
plt.plot(temp80, C80)
plt.ylabel("$C_v$")
plt.title("Model of Phase Transition \n Number of MC cycles = 1E6")
plt.legend(["L=20", "L=40", "L=60", "L=80"], loc=2)

plt.subplot(212)
plt.plot(temp20, X20)
plt.plot(temp40, X40)
plt.plot(temp60, X60)
plt.plot(temp80, X80)
plt.xlabel("Temperature Tk/J")
plt.ylabel("$\\chi$")
plt.legend(["L=20", "L=40", "L=60", "L=80"], loc=2)
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


