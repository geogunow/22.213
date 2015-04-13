import rksolver as rk
import numpy as np
import matplotlib.pyplot as plt
import gc
from utils import *
from math import *

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size' : 18}
plt.rc('font', **font)

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

# set cross sections for fuel3 (3% enriched + 80)
fuel5 = rk.XSdata()
fuel5.setD( [1.3, 0.5] )
fuel5.setSiga( [0.0098, 0.114*0.2 + 0.118*0.8] )
fuel5.setSigs( [[0, 0.022], [0, 0]])
fuel5.setNuSigf( [0.006, 0.195] )
fuel5.setChi( [1.0, 0.0] )

# eigenvalue solver parameters
outer_tol = 10**-6
inner_tol = 10**-5
maxiters = 10000
inner_solver = 1 # optimal SOR

# kinetic parameters
'''
betas = [0.001, 0.002]
halflife = [100, 20]
'''
betas = [0.000218, 0.001023, 0.000605, 0.00131, 0.00220, 0.00060, 0.000540,
        0.000152]
halflife = [55.6, 24.5, 16.3, 5.21, 2.37, 1.04, 0.424, 0.195]
lambdas = log(2) / np.array(halflife)
energy = np.array([10**3, 0.1])
v = 2200 * 100 * np.sqrt(energy / 0.0253)
beta = sum(betas)
chi_d = [1,0]
print 'beta-eff = ', beta

# set kinetic parameters
rkParams = rk.RKdata()
rkParams.setChiD(chi_d)
rkParams.setV(v)
rkParams.setBetas(betas)
rkParams.setLambdas(lambdas)

# describe geometry
mats = []
mats.append([water, fuel1, fuel1, fuel1, water])
mats.append([water, fuel1, fuel1, fuel1, water])
widths = [30.0, 170.0, 20.0, 170.0, 30.0]
nodes = [30, 170, 20, 170, 30]
N = sum(nodes)

# form all mesh
mesh = []
for i in range(2):
    temp = rk.Mesh()
    temp.setMesh(nodes)
    temp.setWidths(widths)
    temp.setMaterials(mats[i])
    temp.setBoundaryConditions(2,2)
    mesh.append(temp)

# set transient
trans = rk.Transient()
trans_t = [0, 2, 4, 10, 12, 50]
trans_m = [1, 1, 0, 0, 1, 1]
trans.setInterpTimes(trans_t)
meshVector = []
for i in xrange(len(trans_t)):
    meshVector.append(mesh[ trans_m[i] ])
trans.setMeshVector(meshVector)

# set time steps
t = np.linspace(0, 50, 100)
trans.setCalcTimes(t)

# solve reactor kinetics problem
result = rk.solveTransientFT(trans, rkParams)
P = np.array(result.getPower())
P = P/P[0]

plt.plot(t, P, 'kx-')
plt.show()
