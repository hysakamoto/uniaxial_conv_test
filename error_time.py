from subprocess import call
from numpy import *
import matplotlib.pyplot as plt

from numerical_integration import *    
from plot_errors import *
from IO_sim import *


n_sim = 10

## load the initial x,y,z values for integration
print "loading initial x,y,z values"
nodes_set = []
elems_set = []
for i in range(1,n_sim):
    print i
    filename = 'time_test/test_{0:d}.feb'.format(i)
    nodes, elems = read_nodes_elems(filename )
    nodes_set.append(nodes)
    elems_set.append(elems)

n_rec_nodes = len(nodes_set[0])
## load the x,y,z,p values for each node
xyzp_set = []
print "loading x,y,z,p values for each node..."
for i in range(1,n_sim):
    print i
    xyzp_set.append( read_xyzp('time_test/test_{0:d}.dat'.format(i), n_rec_nodes, -1 ))

print "space integrate the errors"
Exs = []
Eys = []
Ezs = []
Eps = []
dt = []
for i in range(1,n_sim):
    print "\n", i
    Ex, Ey, Ez, Ep = \
        space_integrate_error( xyzp_set[-1][-1], xyzp_set[i-1][-1], \
                               nodes_set[-1], elems_set[-1]  )
    
    dt.append(0.5/2.0**float(i-1))

        # point_error( xyzp_set[n_sim-1][-1], xyzp_set[i][-1], \
        #                        nodes_set[-1], elems_set[-1]  )

    Exs.append(Ex)
    Eys.append(Ey)
    Ezs.append(Ez)
    Eps.append(Ep)

print "Ex = ", Exs
print "Ey = ", Eys
print "Ex = ", Ezs
print "Ep = ", Eps

plot_errors(Exs, Eys, Ezs, Eps, dt)


# xmin = 0.0
# xmax = 0.5
# ymin = 0.0
# ymax = 0.5
# zmin = 0.0
# zmax = 1.0
# xlen = 0.5
# ylen = 0.5
# zlen = 1.0

# tree = []
# max_depth = 4
# # for i in range(1):
# tree = build_octree(0, tree, elems, nodes, xmin, xmax, ymin, ymax, zmin, zmax, xlen, ylen, zlen, max_depth)


# xyz_range = [[0,0.5], [0,0.5], [0,1.0]]
# sim_num = 2
# max_depth = 5

# tree = space_integrate_error_2 ( xyzp_set[-1][-1], xyzp_set[sim_num][-1], \
#                           nodes_set[-1], nodes_set[sim_num], \
#                           elems_set[-1], elems_set[sim_num], \
#                           max_depth,\
#                           xyz_range )
