from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
from numerical_integration import * 

def plot_errors (Ex, Ey, Ez, Ep, ms):
    "Erros in x, y, z, p (Ex, Ey, Ez, Ep) and mesh size (ms)"

    n_sim = len(Ep)

    # t = [2.0/(2**(i+1)) for i in range(n_sim)]

    ## x-displacement error
    plt.plot(ms,Ex,'-o')
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.title('x-displacement error')
    plt.show()

    ## y-displacement error
    plt.plot(ms,Ey,'-o')
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.title('y-displacement error')
    plt.show()

    ## z-displacement error
    plt.plot(ms,Ez,'-o')
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.title('z-displacement error')
    plt.show()

    ## p-displacement error
    plt.plot(ms,Ep,'-o')
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.title('p error')
    plt.show()


    print "x error convergence = %e" %((log(Ex[-2])-log(Ex[1])) / (log(ms[-2]) - log(ms[1])))
    print "y error convergence = %e" %((log(Ey[-2])-log(Ey[1])) / (log(ms[-2]) - log(ms[1])))
    print "z error convergence = %e" %((log(Ez[-2])-log(Ez[1])) / (log(ms[-2]) - log(ms[1])))
    print "p error convergence = %e" %((log(Ep[-2])-log(Ep[1])) / (log(ms[-2]) - log(ms[1])))
