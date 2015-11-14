import numpy as np
import matplotlib.pyplot as plt

f = "c_ordered_L20t1mcE4.txt"
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
chi_abs    = [] # with <|M|>
a_count    = []
Evariance  = []
E_count = []

for line in data_file:
    line = line.split()
    if len(line)>5:
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
        a_count.append(float(line[10]))
        Evariance.append(float(line[11]))
        if len(line)>12:
            count = []
            for number in line[12:]:
                count.append(float(number))
            E_count.append(count) # nested list
            
def plot_relative_error(E_exp, M_exp, Cv, ch_abs, number_mcs):
    """
    L=2, T=1
    """
    E_comp = np.array(E_exp)
    M_abs_comp = np.array(M_exp)
    Cv_comp = np.array(Cv)
    X_abs_comp = np.array(chi_abs)
    
    # Analytical
    from numpy import exp, cosh, sinh
    N = 4.0;
    temp = 1.0
    E     = -8*sinh(8/temp)/(3+cosh(8/temp))/N
    exp_M     = 0
    Mabs = 2*(exp(8/temp) + 2)/(cosh(8/temp) + 3)/N 
    C_v       = (64./(temp*temp))*(cosh(8/temp)/(3+cosh(8/temp)) - (sinh(8/temp)/(3+cosh(8/temp)))**2)/N
    #Cv      = 64/(temp**2) * (3*cosh(8/temp) + 1)/((cosh(8/temp)+3)**2) #same, differently written
    X         = (8./temp)*(exp(8./temp) + 1)/(3+cosh(8./temp))/N
    Xabs     = (4./(temp*(cosh(8./temp)+3)))*(2*(exp(8./temp) +1 ) - ((exp(8./temp)+2)**2)/(cosh(8./temp)+3))/N

    erE = abs((E_comp-E)/E)
    erM = abs((M_abs_comp-Mabs)/Mabs)
    erCv = abs((Cv_comp-C_v)/C_v)
    erX = abs((X_abs_comp-Xabs)/Xabs)
    plt.figure(1)
    plt.subplot(211)
    plt.plot(number_mcs, erE)
    plt.plot(number_mcs, erM)
    plt.title("Relative error computed vs. analytical \nL=2 T=1.0")
    plt.ylabel("Relative error")
    plt.legend(["$\\langle E\\rangle$", "$\\langle |M| \\rangle$"])

    plt.subplot(212)
    plt.plot(number_mcs, erCv)
    plt.plot(number_mcs, erX)
    plt.ylabel("Relative error")
    plt.xlabel("Time [number of MC cycles]")
    plt.legend(["$Cv$", "$\\chi$"])
    plt.show()

def c_accept():
    number_mcs = np.array(number_mcs)
    count1 = np.array(a_count1)/number_mcs
    count2 = np.array(a_count2)/number_mcs
    
    plt.figure(1)
    plt.subplot(211)
    plt.plot(number_mcs[0:4000], count2[0:4000], 'r-')
    plt.title("Number of accepted spin flips\nL=20")
    plt.legend(["T=2.4\ninit. config.=random"], loc=1)
    plt.grid()
    
    plt.subplot(212)
    plt.plot(number_mcs[0:4000], count1[0:4000], 'g-')
    plt.ylabel("Number of accepted moves")
    plt.legend(["T=1.0\ninit. config.=ordered"], loc=1)
    plt.grid()
    
    plt.xlabel("Time [number of MC cycles]")
    plt.show()
    

def c_exp():
    # Mean energy
    plt.figure(1)
    plt.subplot(211)
    plt.plot(number_mcs, E_exp, '-')
    plt.ylabel("Energy/J")
    plt.title("Expectation values approaching equilibrium\n L=20, T=1.0, init. config. = ordered")
    plt.legend(["$\\langle E\\rangle$"], loc=2)
    plt.grid()
    
    
    # Mean magnetization
    plt.subplot(212)
    plt.plot(number_mcs, M_exp, '-')
    plt.plot(number_mcs, M_abs_exp, '-')
    plt.ylabel("Magnetization")
    plt.legend(["$|\\langle M \\rangle|$","$\\langle |M| \\rangle$"], loc=7)
    plt.grid()
    
    plt.xlabel("Time [number of MC cycles]")
    plt.show()

def d_count():
                    
    n_mcs = number_mcs[0]
    L = L[0]
    energies = np.linspace(-2*L*L, 2*L*L, len(Et1))
    
    Et1 = np.array(Et1)/n_mcs
    Et24 = np.array(Et24)/n_mcs
    print sum(Et1)
    print sum(Et24)
    
    print Evariance
    
    plt.figure(1)
    plt.subplot(211)
    plt.plot(energies, Et24)
    plt.ylabel("P(E)")
    plt.legend(["T=2.4"])
    plt.title("P(E) at L=20")
    
    plt.subplot(212)
    plt.plot(energies, Et1, 'r')
    plt.axis([-810, 800, 0, 0.9])
    plt.legend(["T=1.0"])
    plt.xlabel("E")
    plt.show()
