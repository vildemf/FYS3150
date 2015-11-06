from numpy import *
from matplotlib import pyplot as plt

N = 4.0;
temp = 1.0

#exp_E2   = 16*cosh(8/temp)/(3+cosh(8/temp))
exp_E     = -8*sinh(8/temp)/(3+cosh(8/temp))
exp_M     = 0
exp_M_abs = 2*(exp(8/temp) + 2)/(cosh(8/temp) + 3) 
C_v       = (64./(temp*temp))*(cosh(8/temp)/(3+cosh(8/temp)) - (sinh(8/temp)/(3+cosh(8/temp)))**2)
#C_v      = 64/(temp**2) * (3*cosh(8/temp) + 1)/((cosh(8/temp)+3)**2) #same, differently written
X         = (8./temp)*(exp(8./temp) + 1)/(3+cosh(8./temp))
X_abs     = (4./(temp*(cosh(8./temp)+3)))*(2*(exp(8./temp) +1 ) - ((exp(8./temp)+2)**2)/(cosh(8./temp)+3))



#print exp_E2/N
print exp_E/N
print exp_M/N
print exp_M_abs/N
print C_v/N
print X/N
print X_abs/N


temp = linspace(-5, 5, 100)
C_v = (64/(temp*temp))*(cosh(8*temp)/(3+cosh(8*temp)) - (sinh(8*temp)/(3+cosh(8*temp)))**2)
#plt.plot(temp, C_v)
#plt.show()
