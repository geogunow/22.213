import rksolver as rk
import numpy as np
import matplotlib.pyplot as plt
from utils import *
from math import *

# number of groups
G = 2

# set cross sections for water
water = rk.XSdata()
water.setD( [1.5, 0.5] )
water.setSiga( [0.0002, 0.010] )
water.setSigs( [[0, 0.032], [0, 0]] )
water.setNuSigf( [0, 0] )
water.setChi( [1, 0] )

# set cross sections for fuel1 (3% enriched)
fuel1 = rk.XSdata()
fuel1.setD( [1.3, 0.5] )
fuel1.setSiga( [0.0098, 0.114] )
fuel1.setSigs( [[0, 0.022], [0, 0]])
fuel1.setNuSigf( [0.006, 0.195] )
fuel1.setChi( [1.0, 0.0] )

# set cross sections for fuel2 (4% enriched)
fuel2 = rk.XSdata()
fuel2.setD( [1.3, 0.5] )
fuel2.setSiga( [0.0098, 0.114] )
fuel2.setSigs( [[0, 0.022], [0, 0]])
fuel2.setNuSigf( [0.008, 0.23] )
fuel2.setChi( [1.0, 0.0] )

# set cross sections for black absorber
black = rk.XSdata()
black.setD( [1.3, 0.5] )
black.setSiga( [999.9, 999.9] )
black.setSigs( [[0, 0.020], [0, 0]])
black.setNuSigf( [0, 0] )
black.setChi( [1.0, 0.0] )

# set cross sections for fuel3 (3% enriched + rod1)
fuel3 = rk.XSdata()
fuel3.setD( [1.3, 0.5] )
fuel3.setSiga( [0.0098, 0.118] )
fuel3.setSigs( [[0, 0.022], [0, 0]])
fuel3.setNuSigf( [0.006, 0.195] )
fuel3.setChi( [1.0, 0.0] )

# set cross sections for fuel3 (3% enriched + rod2)
fuel4 = rk.XSdata()
fuel4.setD( [1.3, 0.5] )
fuel4.setSiga( [0.0150, 0.145] )
fuel4.setSigs( [[0, 0.022], [0, 0]])
fuel4.setNuSigf( [0.006, 0.195] )
fuel4.setChi( [1.0, 0.0] )

# describe geometry
mats = [water, fuel1, fuel3, fuel1, water]
widths = [30.0, 170.0, 20.0, 170.0, 30.0]
nodes = [3, 17, 2, 17, 3]
N = sum(nodes)

# eigenvalue solver parameters
outer_tol = 10**-6
inner_tol = 10**-5
maxiters = 10000
inner_solver = 1 # optimal SOR

# kinetic parameters
beta_i = [0.000218, 0.001023, 0.000605, 0.00131, 0.00220, 0.00060, 0.000540,
        0.000152]
halflife = [55.6, 24.5, 16.3, 5.21, 2.37, 1.04, 0.424, 0.195]
lambda_i = log(2) / np.array(halflife)
velocity1 = 2200 * 100 * sqrt(10**3 / 0.0253)
velocity2 = 2200 * 100 * sqrt(0.1 / 0.0253)
beta = sum(beta_i)
print 'beta-eff = ', beta


'''
Parts D and E
'''
# form mesh
mesh = rk.Mesh()
mesh.setMesh(nodes)
mesh.setWidths(widths)
mesh.setMaterials(mats)
mesh.setBoundaryConditions(2,2)

# solve eigenvalue problem
solution = rk.solveDiffusion(mesh, outer_tol, inner_tol, maxiters,
        inner_solver)

