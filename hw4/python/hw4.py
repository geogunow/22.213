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
plt.rc('legend', **{'fontsize':18})

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

# set cross sections for fuel3 (3% enriched + rod1)
fuel3 = rk.XSdata()
fuel3.setD( [1.3, 0.5] )
fuel3.setSiga( [0.0098, 0.118] )
fuel3.setSigs( [[0, 0.022], [0, 0]])
fuel3.setNuSigf( [0.006, 0.195] )
fuel3.setChi( [1.0, 0.0] )

# eigenvalue solver parameters
outer_tol = 10**-6
inner_tol = 10**-5
maxiters = 10000
inner_solver = 1 # optimal SOR

# kinetic parameters
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
mats.append([water, fuel1, fuel3, fuel1, water])

widths = [30.0, 170.0, 20.0, 170.0, 30.0]
nodes = [12, 68, 8, 68, 12]
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

t = np.linspace(0,50,5001)
trans.setCalcTimes(t)

'''
Problem A
'''
mesh_out = mesh[0]
mesh_in = mesh[1]

ans_out = rk.solveCritical(mesh_out, 10**-8, 10**-8, 1000, 2)
ans_in = rk.solveCritical(mesh_in, 10**-8, 10**-8, 1000, 2)

# plot results
x = mesh_out.getX()
plt.plot(x, ans_out.getFlux(1),'bx-')
plt.plot(x, ans_out.getFlux(2),'gx-')
plt.xlabel('Distance (cm)')
plt.ylabel('Normalized Flux')
plt.show()

x = mesh_in.getX()
plt.plot(x, ans_in.getFlux(1),'bx-')
plt.plot(x, ans_in.getFlux(2),'gx-')
plt.xlabel('Distance (cm)')
plt.ylabel('Normalized Flux')
plt.show()


'''
Problem B
'''
result = rk.solveTransient(trans, rkParams)
P_ref = np.array(result.getPower())
P_ref /= P_ref[0]
plt.plot(t,P_ref,'k.-')
plt.xlabel('Time (s)')
plt.ylabel('Normalized Power')
plt.show()

'''
Problem C
'''
# part 1
result = rk.solvePKE(trans, rkParams, 0, 0, 0)
P = np.array(result.getPower())
plt.plot(t,P_ref,'k.-')
plt.plot(t,P/P[0],'r.-')
result = rk.solvePKESimple(trans, rkParams, 0, 0, 0)
P = np.array(result.getPower())
plt.plot(t,P/P[0],'g.-')
plt.legend(['Reference', 'Point Kinetics Two Group', 'Point Kinetics One Group'])
plt.xlabel('Time (s)')
plt.ylabel('Normalized Power')
plt.show()

# part 2
result = rk.solvePKE(trans, rkParams, 0, 1, 0)
P = np.array(result.getPower())
plt.plot(t,P_ref,'k.-')
plt.plot(t,P/P[0],'r.-')
result = rk.solvePKESimple(trans, rkParams, 0, 1, 0)
P = np.array(result.getPower())
plt.plot(t,P/P[0],'g.-')
plt.legend(['Reference', 'Point Kinetics Two Group', 'Point Kinetics One Group'],
        prop={'size':18})
plt.xlabel('Time (s)')
plt.ylabel('Normalized Power')
plt.show()

# part 3
result = rk.solvePKE(trans, rkParams, 3, 0, 0)
P = np.array(result.getPower())
plt.plot(t,P_ref,'k.-')
plt.plot(t,P/P[0],'r.-')
result = rk.solvePKESimple(trans, rkParams, 3, 0, 0)
P = np.array(result.getPower())
plt.plot(t,P/P[0],'g.-')
plt.legend(['Reference', 'Point Kinetics Two Group', 'Point Kinetics One Group'],
        prop={'size':18})
plt.xlabel('Time (s)')
plt.ylabel('Normalized Power')
plt.show()

# part 4
result = rk.solvePKE(trans, rkParams, 3, 1, 3)
P = np.array(result.getPower())
plt.plot(t,P_ref,'k.-')
plt.plot(t,P/P[0],'r.-')
result = rk.solvePKESimple(trans, rkParams, 3, 1, 3)
P = np.array(result.getPower())
plt.plot(t,P/P[0],'g.-')
plt.legend(['Reference', 'Point Kinetics Two Group', 'Point Kinetics One Group'],
        prop={'size':18})
plt.xlabel('Time (s)')
plt.ylabel('Normalized Power')
plt.show()

# part 5
result = rk.solvePKE(trans, rkParams, 0, 1, 3)
P = np.array(result.getPower())
plt.plot(t,P_ref,'k.-')
plt.plot(t,P/P[0],'r.-')
result = rk.solvePKESimple(trans, rkParams, 0, 1, 3)
P = np.array(result.getPower())
plt.plot(t,P/P[0],'g.-')
plt.legend(['Reference', 'Point Kinetics Two Group', 'Point Kinetics One Group'],
        prop={'size':18})

plt.xlabel('Time (s)')
plt.ylabel('Normalized Power')
plt.show()

# part 6
result = rk.solvePKE(trans, rkParams, 3, 1, 0)
P = np.array(result.getPower())
plt.plot(t,P_ref,'k.-')
plt.plot(t,P/P[0],'r.-')
P = np.array(result.getPower())
plt.plot(t,P/P[0],'g.-')
plt.legend(['Reference', 'Point Kinetics Two Group', 'Point Kinetics One Group'],
        prop={'size':18})
plt.xlabel('Time (s)')
plt.ylabel('Normalized Power')
plt.show()
