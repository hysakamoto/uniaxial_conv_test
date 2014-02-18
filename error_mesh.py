from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
from numerical_integration import *    
from plot_errors import *
from IO_sim import *


n_sim = 5
xyz_range = [0.0, 0.5, 0.0, 0.5, 0.0, 1.0]
xlen = xyz_range[1] - xyz_range[0]
ylen = xyz_range[3] - xyz_range[2]
zlen = xyz_range[5] - xyz_range[4]

## load the initial x,y,z values for integration
print "loading initial x,y,z values"
nodes_set = []
elems_set = []
for i in range(n_sim):
    print i
    filename = 'mesh_test/test_{0:d}.feb'.format(i)
    nodes, elems = read_nodes_elems(filename )
    nodes_set.append(nodes)
    elems_set.append(elems)

## load the x,y,z,p values for each node
xyzp_set = []
end_step = 102
print "loading x,y,z,p values for each node..."
for i in range(n_sim):
    print i
    n_rec_nodes = len(nodes_set[i]) 
    # n_rec_nodes = 8
    xyzp_set.append( read_xyzp('mesh_test/test_{0:d}.dat'.format(i), n_rec_nodes, end_step ))


print "space integrate the errors..."
Exs = []
Eys = []
Ezs = []
Eps = []
for i in range(0,n_sim):
    print "sim: ", i
    max_depth = int(float(i)/2.0)
    Ex, Ey, Ez, Ep = \
        space_integrate_error_2( xyzp_set[n_sim-1][-1], xyzp_set[i][-1], 
                                 nodes_set[-1], nodes_set[i], \
                                 elems_set[-1], elems_set[i],
                                 max_depth, xyz_range )
    # Ex, Ey, Ez, Ep = \
    #     point_error( xyzp_set[n_sim-1][-1], xyzp_set[i][-1], \
    #                            nodes_set[-1], elems_set[-1]  )

    print "\nerror: " %i, Ex, Ey, Ez, Ep

    Exs.append(Ex)
    Eys.append(Ey)
    Ezs.append(Ez)
    Eps.append(Ep)


ms = [(float(j)/2.0+1.0)**2.0 for j in range(n_sim)]
plot_errors (Exs[0:-1], Eys[0:-1], Ezs[0:-1], Eps[0:-1], ms[0:-1])
