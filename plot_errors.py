from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
from numerical_integration import * 

def plot_errors (Ex, Ey, Ez, Ep, ms):
    "Erros in x, y, z, p (Ex, Ey, Ez, Ep) and mesh size (ms)"

    # Ex = [0.020414517073643605, 0.020330688276197547, 0.020323384507120797, 0.02012798103791228, 0.01976711641518359, 0.017799423297099954, 0.016408267182855197, 0.01608153814148484, 0.0] 

    # Ey = [0.020414517073973647, 0.020330688276524414, 0.02032338450744697, 0.020127981038238443, 0.01976711641550977, 0.01779942329742615, 0.01640826718318138, 0.01608153814181077, 0.0] 

    # Ez = [0.06866535002118805, 0.06855187176190802, 0.06850911923187693, 0.06820883845892774, 0.06765841582990413, 0.06439373524329027, 0.0611560847697454, 0.05935122165150506, 0.0]

    # Ep = [0.0163707851121223, 0.016421242916140777, 0.016593684021803858, 0.016447993908707662, 0.016335161251599624, 0.015193080261123704, 0.01247765424572639, 0.008164325930323485, 0.0]


    # Ex = [1.2771335333275344e-07, 7.3918602518839458e-08, 3.9923598328698869e-08, 2.0112131963177344e-08, 9.8210222883138284e-09, 4.5316633568016091e-09, 1.9289603079417727e-09, 5.9731793810663506e-10, 0.0]

    # Ey = [1.2771211377949971e-07, 7.3917362991525031e-08, 3.9922358792818663e-08, 2.0110892428757469e-08, 9.8197827523769786e-09, 4.5304238207506455e-09, 1.9277207705268495e-09, 5.9607840079266524e-10, 0.0] 

    # Ez = [8.135300523770315e-07, 4.7193165595673021e-07, 2.5241007866780187e-07, 1.2754969567606732e-07, 6.2419290243039261e-08, 2.8949667806643863e-08, 1.2481074235007374e-08, 4.0448076519544283e-09, 0.0] 

    # Ep = [0.064870094336657066, 0.037215213767036244, 0.020218030860998291, 0.010160672988531435, 0.0049576840192060868, 0.0022835367990901223, 0.00096769615851211314, 0.0002949048956693337, 0.0]

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
