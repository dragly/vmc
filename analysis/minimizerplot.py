from enthought.mayavi.mlab import *
from numpy import *

param0 = loadtxt("../vmc-debug-mpi/parameters0.dat")
param1 = loadtxt("../vmc-debug-mpi/parameters1.dat")
data = loadtxt("../vmc-debug-mpi/energies.dat")

surf(param0, param1, data, warp_scale='auto')

axes()
