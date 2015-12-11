from numpy import *
from matplotlib.pyplot import *
from os import system
import sys


cpp_program = "build-diffusion-Desktop_Qt_5_4_1_GCC_64bit-Debug/diffusion"
T_list =        ['0.001', '0.01', '0.08', '0.3']#sys.argv[2].split(',')
d =             1 #float(sys.argv[3])
Nt =            50000  #int(sys.argv[4])
Nx =            200 #int(sys.argv[5])
forward =       0#int(sys.argv[6])
backward =      0#int(sys.argv[7])
crank =         0#int(sys.argv[8])
store_x_not_t = 1#int(sys.argv[9])

if forward==1:
    method='forward'
    titlename='Explicit Forward Euler'
elif backward==1:
    method='backward'
    titlename='Implicit Backward Euler'
elif crank==1:
    method='crank'
    titlename='Implicit Crank-Nicolson'


# "./build-diffusion-Desktop_Qt_5_4_1_GCC_64bit-Debug/diffusion argtest.txt 0.5 1 50000 200 1 0 0 1"

    
    

def u_analytical(x, t):
    v = 0
    for n in range(1,50):
        v += (1./n)*sin(n*pi*x)*exp(-t*(n*pi)**2)
    v *= -2/pi
    u = array(v) + 1 - x
    return u



def plot_deltaxerror(T):
    T1 = T.split('.')[0]
    T2 = T.split('.')[1]
    T = float(T)
    Nx_list = linspace(120, 280, 10)
    forward_points = []
    backward_points = []
    crank_points = []
    dx_list = []
    for Nx in Nx_list:
        x = linspace(0, d, Nx+1)
        dx_list.append(x[1]-x[0])
        u_an = u_analytical(x, T)[1]
        for i in range(3):
            if i==1:
                method='forward'
                titlename='Explicit Forward Euler'
                forward = 1
                backward = 0
                crank = 0
                filename = "xerror_%s_T%sp%sNt%dNx%d.txt" % (method, T1, T2, Nt, Nx)         
                system("./%s %s %f %f %d %d %d %d %d %d" % \
                       (cpp_program, filename, T, d, Nt, Nx, forward, backward, crank, store_x_not_t))
        
                datafile = open(filename)
                forward_points.append(abs(float(datafile.readlines()[1])+1-x[1] - u_an))
            elif i==2:
                method='backward'
                titlename='Implicit Backward Euler'
                backward = 1
                forward = 0
                crank = 0
                filename = "xerror_%s_T%sp%sNt%dNx%d.txt" % (method, T1, T2, Nt, Nx)         
                system("./%s %s %f %f %d %d %d %d %d %d" % \
                       (cpp_program, filename, T, d, Nt, Nx, forward, backward, crank, store_x_not_t))
        
                datafile = open(filename)
                backward_points.append(abs(float(datafile.readlines()[1])+1-x[1] - u_an))
            else: # i==3:
                method='crank'
                titlename='Implicit Crank-Nicolson'
                crank = 1
                backward = 0
                forward = 0
                filename = "xerror_%s_T%sp%sNt%dNx%d.txt" % (method, T1, T2, Nt, Nx)         
                system("./%s %s %f %f %d %d %d %d %d %d" % \
                       (cpp_program, filename, T, d, Nt, Nx, forward, backward, crank, store_x_not_t))
        
                datafile = open(filename)
                crank_points.append(abs(float(datafile.readlines()[1])+1-x[1] - u_an))
    #flist = array(forward_points)
    #blist = array(
    semilogy(list(reversed(dx_list)), list(reversed(forward_points)))
    semilogy(list(reversed(dx_list)), list(reversed(backward_points)))
    semilogy(list(reversed(dx_list)), list(reversed(crank_points)))
    legend(['Forward', 'Backward', 'Crank'])
    show()
#plot_deltaxerror(T_list[0])


def plot_xerror(T, abs_error):
    x = linspace(0, d, Nx+1)
    dx =  x[1]-x[0]
    T1 = T.split('.')[0]
    T2 = T.split('.')[1]
    T=float(T)
    dt = linspace(0, T, Nt+1)[1] - linspace(0, T, Nt+1)[0]
    legends= []
    for i in range(3):
        if i==1:
            method='forward'
            titlename='Explicit Forward Euler'
            forward = 1
            backward = 0
            crank = 0
        elif i==2:
            method='backward'
            titlename='Implicit Backward Euler'
            backward = 1
            forward = 0
            crank = 0
        else: # i==3:
            method='crank'
            titlename='Implicit Crank-Nicolson'
            crank = 1
            backward = 0
            forward = 0
        filename = "xerror_%s_T%sp%sNt%dNx%d.txt" % (method, T1, T2, Nt, Nx)         
        system("./%s %s %f %f %d %d %d %d %d %d" % \
               (cpp_program, filename, T, d, Nt, Nx, forward, backward, crank, store_x_not_t))
        
        datafile = open(filename)
        v = []
        for line in datafile:
            v.append(float(line))
        u_num = array(v) + 1 - x
        u_an = u_analytical(x, T)
        if abs_error:
            error = abs(u_num-u_an)
            error_name = "Absolute error"
        else:
            error = abs(u_num-u_an)/u_an
            error_name = "Relative error"
        plot(x, error*1E4)
        legends.append("%s" % (titlename))
    override = {
    'fontsize'            : 'large',
    'verticalalignment'   : 'baseline',
    'horizontalalignment' : 'center'
    }
    xlabel("x", override)
    ylabel("%s $\\times 1E-4$" % (error_name), override)
    title("%s (steady state T=%.1f) \n Nt=%d $\\Delta t$=%.0e Nx=%d $\\Delta x$=%.0e" % (error_name, T,Nt,dt,Nx,dx))
    legend(legends, loc=6)
    show()

#plot_xerror(T_list[0], True)

def plot_solution(T_list, d, Nt, Nx, method, cpp_program, titlename):
    legends = []
    x = linspace(0, d, Nx+1)
    dx = x[1]-x[0]
    for T in T_list:
        T1 = T.split('.')[0]
        T2 = T.split('.')[1]
        T=float(T)
        dt = linspace(0, T, Nt+1)[1] - linspace(0, T, Nt+1)[0]
        filename = "%s_T%sp%sNt%dNx%d.txt" % (method, T1, T2, Nt, Nx) 
        
        #system("./%s %s %f %f %d %d %d %d %d %d" % \
        #   (cpp_program, filename, T, d, Nt, Nx, forward, backward, crank, store_x_not_t))

        datafile = open(filename)
        v = []
        for line in datafile:
            v.append(float(line))        
        u = array(v) + 1 - x
        plot(x, u)
        legends.append('T=%.3f $\\Delta t$=%.0e' % (T, dt))
    #override = {
    #'fontsize'            : 'large',
    #'verticalalignment'   : 'baseline',
    #'horizontalalignment' : 'center'
    #}
    xlabel("x", override)
    ylabel("u(x, t=T)", override)
    title("%s \n Nt=%d Nx=%d $\\Delta x$=%.0e" % (titlename, Nt, Nx, dx))
    legend(legends)
    #show()


titlename = 'h'

plot_solution(T_list, d, Nt, Nx, forward, cpp_program, titlename)    
plot_solution(T_list, d, Nt, Nx, backward, cpp_program, titlename) 
plot_solution(T_list, d, Nt, Nx, crank, cpp_program, titlename) 

override = {
    'fontsize'            : 'large',
    'verticalalignment'   : 'baseline',
    'horizontalalignment' : 'center'
}
xlabel("x", override)
ylabel("u(x, t=T)", override)
title("%s \n Nt=%d Nx=%d $\\Delta x$=%.0e" % (titlename, Nt, Nx, dx))
legend(legends)
show()
