uniaxial_conv_test
==================
The convergence analysis of the finite element simulation of poroviscleastic problem for a uniaxial tension case

### Prerequisites
* [FEBio](http://www.febio.org/) --> for running finite element simulation
* [gmsh](http://geuz.org/gmsh/) --> for generating mesh in .msh format.
* lxml --> for creating FEBio input file (xml format)
* numpy, matplotlib

### Simulation file generation
* make_test.py --> generates a test simulation file.
* make_mesh_test.py --> generates mesh convergence test simulation files.
* make_time_test.py --> generates time convergence test simulation files.

### Run simulations using FEBio
Something like
```bash
febio test_mesh/test_0.feb
```

### Convergence analysis
* error_mesh.py --> convergence analysis on space
* error_time.py --> convergence analysis on time
