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

# p =  [0.2, 0.121, 0.923]

# sim_num = 1
# tree = []
# for i in range(len(elems_set[sim_num])):
#     max_depth = sim_num
#     tree = build_octree(i, tree, elems_set[sim_num], nodes_set[sim_num], xyz_range, \
#                         xlen, ylen, zlen, max_depth)


# # inv_elem_n, xis = \
# #         inverse_physical_to_master( p, xyz_range, tree, \
# #                                     nodes_set[sim_num], elems_set[sim_num], \
# #                                     sim_num, max_depth )

# # xv, yv, zv = get_xyzp_values( nodes_set[sim_num], elems_set[sim_num], \
# #                               xyzp_set[sim_num][-1], inv_elem_n, xis )

# # print xv, yv, zv

# E = \
#     [ [3.93978590604e-06, 3.93972449486e-06, 7.99170400625e-06, 2.05815737354], \
#       [3.16419307349e-06, 3.16415434278e-06, 4.98769053813e-06, 0.875692828492], \
#       [2.22100313638e-06, 2.22097286242e-06, 3.08040418243e-06, 0.238793716379], \
#       [1.62122257528e-06, 1.62119918274e-06, 2.2181971831e-06, 0.109050169694], \
#       [1.01783956895e-06, 1.01782675289e-06, 1.37422692692e-06, 0.0485364031153], \
#       [6.47909994935e-07, 6.47906009e-07, 8.60987329556e-07, 0.0281059042658], \
#       [3.61366011932e-07, 3.61370968771e-07, 4.67053907536e-07, 0.0168283649123], \
#       [1.96133379742e-07, 1.96142671432e-07, 2.46811745349e-07, 0.01221449734] ]

# E = array(E).T

# ms = [1.0/(float(j)/2.0+1.0)**2.0 for j in range(len(E[0]))]
# plot_errors (E[:][0], E[:][1], E[:][2], E[:][3], ms)
