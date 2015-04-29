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

# desribe transients
trans1 = rk.Transient()
trans2 = rk.Transient()
for trans in [trans1, trans2]:

    # define transient specifice parameters such as rodded material
    if trans == trans1:
        rod = fuel3
        trans_t = [0, 2, 4, 10, 12, 50]
    else:
        rod = fuel4
        trans_t = [0, 0.1, 0.2, 0.5, 1.5, 2.0]

    # describe geometry
    mats = []
    mats.append([water, fuel1, fuel1, fuel1, water])
    mats.append([water, fuel1, rod, fuel1, water])
    widths = [30.0, 170.0, 20.0, 170.0, 30.0]
    nodes = [3, 17, 2, 17, 3]
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

    # define transient
    trans_m = [1, 1, 0, 0, 1, 1]
    trans.setInterpTimes(trans_t)
    meshVector = []
    for i in xrange(len(trans_t)):
        meshVector.append(mesh[ trans_m[i] ])
    trans.setMeshVector(meshVector)

'''
Part A
'''
'''
# compute all cases for global FT
time_scales = [50, 2]
transients = [trans1, trans2]
for tmax, trans in zip(time_scales, transients):

    # define tolerance and timestep criteria
    tolerance_criteria = [1e-2, 1e-4, 1e-6]
    timeStep_criteria = [1e-1, 1e-2, 1e-3, 1e-4] 

    # loop over convergence/timestep criteria
    for tolerance in tolerance_criteria:
        
        # initialize figures
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        leg = []

        for timeStep in timeStep_criteria:
            
            # setup transient
            t = np.linspace(0, tmax, int(tmax/timeStep) + 1)
            trans.setCalcTimes(t)
            trans.setTolerance(tolerance)
            
            # solve transient
            result = rk.solveTransientFT(trans, rkParams, 1)
            P = np.array(result.getPower())
            iters = np.array(result.getInnerSteps())
            
            # plot results
            ax1.semilogy(t, P/P[0])
            ax2.semilogy(t[1:], iters)

            # add legend description
            leg.append("$\Delta t$ = " + str(timeStep))
    
        # label axes
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Normalized Power')
        ax1.legend(leg)
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Iteration Count')
        ax2.legend(leg)

        # display plots
        plt.show()
'''

'''
Part B
'''
# compute all cases for global FT
time_scales = [50, 2]
transients = [trans1, trans2]
for tmax, trans in zip(time_scales, transients):

    # define tolerance and timestep criteria
    tolerance = 1e-6
    timeStep = 1e-2

    # setup transient
    t = np.linspace(0, tmax, int(tmax/timeStep) + 1)
    trans.setCalcTimes(t)
    trans.setTolerance(tolerance)
    
    # try global (1) and local (2) frequencies
    for typ in [1,2]:
        
        # solve transient
        result = rk.solveTransientFT(trans, rkParams, typ)
        P = np.array(result.getPower())

        # plot normalized power
        plt.plot(t, P/P[0])    

    # label plot
    plt.xlabel('Time (s)')
    plt.ylabel('Normalized Power')
    plt.legend(['Global $\omega$','Local $\omega$'])
    plt.show()

