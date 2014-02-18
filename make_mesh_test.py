from write_simulation_mesh import *
from subprocess import call
import os.path
from IO_sim import *

n_sim = 12

id=1;
sim_name = "mesh_test" 
name = "test"

sim_time = 0.2
step_size = 0.5
time_steps = int(sim_time/step_size)

min_dtmax = 0.1

phi0 = 0.5
density = 1.0


c1 = 100
c2 = 0.0
k = c1*1000

ksi = 0.01
beta = 2.0

g1 = 1.0
t1 = 10.0

perm = 0.000001

p = 0
l = 2.2
L = 2.0
pressure_top = -(-p+2*c1*(l/L)*(l/L))
pressure_side = -(-p+2*c1*(L/l))

mu = 2.0*c1
lm = 1.5 # uniaxial tension
pressure_top = mu*(lm-1/(lm*lm))
lm2 = lm**(-1.0/2.0)

# pressure_top = 577.778
pressure_side = 0.0

print "traction = {0:e}".format(pressure_top)
print "lm1 = {0:e}".format(lm)
print "lm2 = {0:e}".format(lm2)

sliding_penalty = 0.0

sim_time = 0.2
step_size = 2.0/(2**(10))
time_steps = int(sim_time/step_size)


#### write the simulation files
id = 1
for i in range(n_sim):

    name = "test_%d" %i
    write_geo_file(i, (float(i)/2.0+1.0)**2.0)
    call(["/Applications/Gmsh.app/Contents/MacOS/gmsh", "-3", "box.geo"])
    call(["cp", "box.msh", "{0:s}/box_{1:d}.msh".format(sim_name, i)])
    call(["cp", "box.geo", "{0:s}/box_{1:d}.geo".format(sim_name, i)])

    if (not os.path.exists("{0:s}/cell_{0:s}_{1:d}.feb".format(sim_name,id))):
		print "writing simulation {0:d}".format(id)
                simulation(id, sim_name, name, pressure_top, pressure_side,  \
                           time_steps, min_dtmax, step_size, \
                           phi0, density, \
                           c1, c2, k, \
                           ksi, beta, \
                           g1, t1, \
                           perm, \
                           sliding_penalty\
                )
    id = id+1

#### Run the simulation
# FEBio = "/h2/ysakamoto/Desktop/ma_febio/FEBio/febio.lnx64"
# FEBio = "./febio"

# for i in range(1,id):
# 	print "running simulation {0:d}".format(i)
# 	# if file not exist
# 	if (not os.path.isfile("{0:s}/cell_{0:s}_{1:d}.xplt".format(sim_name,i))):
# 		call([FEBio, "{0:s}/cell_{0:s}_{1:d}.feb".format(sim_name,i)])
