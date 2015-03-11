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

'''
Problem A
'''
# describe geometry
mats = []
mats.append([water, fuel1, fuel1, fuel1, water])
mats.append([water, fuel1, fuel3, fuel1, water])
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

# calculate static rod worth
sol1 = rk.solveCritical(mesh[1], 10**-6, 10**-6, 10**6, 2)
sol2 = rk.solveCritical(mesh[0], 10**-6, 10**-6, 10**6, 2)

print "Rod worth part B = ", (sol2.keff/sol1.keff - 1)/beta

# set transient
trans = rk.Transient()
trans_t = [0, 2, 4, 10, 12, 50]
trans_m = [1, 1, 0, 0, 1, 1]
trans.setInterpTimes(trans_t)
meshVector = []
for i in xrange(len(trans_t)):
    meshVector.append(mesh[ trans_m[i] ])
trans.setMeshVector(meshVector)

'''
# cycle through number of time steps
min_steps = 10
max_steps = 1000
step = 10
tsteps = np.array(range(min_steps, max_steps + step, step))
for i in range(len(tsteps)):
    tsteps[i] = 5 * int(round(tsteps[i]/5)) + 1

peak_power = []
final_power = []
for npts in tsteps:
    
    print "Solving transient with ", npts, " time steps"

    # set time steps
    t = np.linspace(0,50,npts)
    trans.setCalcTimes(t)

    # solve reactor kinetics problem
    result = rk.solveTransient(trans, rkParams)
    P = np.array(result.getPower())
    P = P/P[0]

    # record power peak, final power
    peak_power.append(max(P))
    final_power.append(P[-1])

# plot error in peak power and final power
ref_peak = peak_power[-1]
ref_final = final_power[-1]
peak_power = np.array(peak_power)
final_power = np.array(final_power)
diff_peak = 100 * abs(peak_power - ref_peak) / ref_peak
diff_final = 100 * abs(final_power - ref_final) / ref_final
plt.semilogy(tsteps, diff_peak, 'kx-')
plt.semilogy(tsteps, diff_final, 'rx-')
plt.xlabel('Number of Time Steps')
plt.ylabel('Error in Power (%)')
plt.legend(['Peak Power', 'Final Power'])
plt.show()

# find number of steps to converge to 1% peak power
peak_steps = tsteps[-1]
for i, error in enumerate(diff_peak):
    if error < 1:
        peak_steps = tsteps[i]
        break;
print "peak power convergence delta-t = ", 50.0/(peak_steps-1)

# find number of steps to converge to 1% final power
for i, error in enumerate(diff_final):
    if error < 1:
        final_steps = tsteps[i]
        break;
print "final power convergence delta-t = ", 50.0/(final_steps-1)

# calculate convergence for problem requested time steps
powers = range(7)
mult = 2**np.array(powers)
peak_power = []
final_power = []
for i, m in enumerate(mult):

    npts = (peak_steps-1)/m + 1
    print "Solving transient with ", npts, " time steps"
    
    # set time steps
    t = np.linspace(0,50,npts)
    trans.setCalcTimes(t)

    # solve reactor kinetics problem
    result = rk.solveTransient(trans, rkParams)
    P = np.array(result.getPower())
    P = P/P[0]

    # record power peak
    peak_power.append(max(P))
    final_power.append(P[-1])

    # plot converged core power vs time
    if i == 0:
        plt.plot(t, P/P[0], 'k.-')
        plt.xlabel('Time (s)')
        plt.ylabel('Relative Power')
        plt.show()

# plot error in peak power and final power
tsteps = 50.0 / ( (peak_steps-1) / np.array(mult) )
peak_power = np.array(peak_power)
final_power = np.array(final_power)
diff_peak = 100 * abs(peak_power - ref_peak) / ref_peak
diff_final = 100 * abs(final_power - ref_final) / ref_final
plt.semilogy(tsteps, diff_peak, 'kx-')
plt.semilogy(tsteps, diff_final, 'rx-')
plt.xlabel('Time Step Size (s)')
plt.ylabel('Error in Power (%)')
plt.legend(['Peak Power', 'Final Power'], loc=5)
plt.show()
'''
peak_steps = 270
# evaluate fractional error in peak, final powers vs spatial mesh
spacing = [1,2,5,10]

peak_power = []
final_power = []
for dx in spacing:

    print "Solving transient for ", dx, " cm mesh spacing"
    
    # form all new mesh
    npts = [node/dx for node in nodes]
    mesh = []
    for i in range(2):
        temp = rk.Mesh()
        temp.setMesh(npts)
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
    t = np.linspace(0, 50, peak_steps)
    trans.setCalcTimes(t)

    # solve reactor kinetics problem
    result = rk.solveTransient(trans, rkParams)
    P = np.array(result.getPower())
    P = P/P[0]

    # record power peak, final power
    peak_power.append(max(P))
    final_power.append(P[-1])

# plot errors in peak, final power
ref_peak = peak_power[0]
ref_final = final_power[0]
peak_power = np.array(peak_power)
final_power = np.array(final_power)
diff_peak = 100 * abs(peak_power - ref_peak) / ref_peak
diff_final = 100 * abs(final_power - ref_final) / ref_final
plt.plot(spacing, diff_peak, 'kx-')
plt.plot(spacing, diff_final, 'rx-')
plt.xlabel('Mesh Spacing (cm)')
plt.ylabel('Error in Power (%)')
plt.legend(['Peak Power', 'Final Power'])
plt.show()

'''
Problem B
'''
# describe geometry
mats = []
mats.append([water, fuel1, fuel1, fuel1, water])
mats.append([water, fuel1, fuel4, fuel1, water])
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

# calculate static rod worth
sol1 = rk.solveCritical(mesh[1], 10**-6, 10**-6, 10**6, 2)
sol2 = rk.solveCritical(mesh[0], 10**-6, 10**-6, 10**6, 2)
print "Rod worth part B = ", (sol2.keff/sol1.keff - 1)/beta

# set transient
trans = rk.Transient()
trans_t = [0, 0.1, 0.2, 0.5, 1.5, 2.0]
trans_m = [1, 1, 0, 0, 1, 1]
trans.setInterpTimes(trans_t)
meshVector = []
for i in xrange(len(trans_t)):
    meshVector.append(mesh[ trans_m[i] ])
trans.setMeshVector(meshVector)

# cycle through number of time steps
'''
min_steps = 10**4
max_steps = 2*10**5
step = 10**4
tsteps = np.array(range(min_steps, max_steps + step, step))
for i in range(len(tsteps)):
    tsteps[i] = 4 * int(round(tsteps[i]/4)) + 1

peak_power = []
final_power = []
for i, npts in enumerate(tsteps):
    
    print "Solving transient with ", npts, " time steps"

    # set time steps
    t = np.linspace(0,2,npts)
    trans.setCalcTimes(t)

    # solve reactor kinetics problem
    result = rk.solveTransient(trans, rkParams)
    P = np.array(result.getPower())
    P = P/P[0]
    del result
    gc.collect()

    # record power peak, final power
    peak_power.append(max(P))
    final_power.append(P[-1])

# plot error in peak power and final power
ref_peak = peak_power[-1]
ref_final = final_power[-1]
peak_power = np.array(peak_power)
final_power = np.array(final_power)
diff_peak = 100 * abs(peak_power - ref_peak) / ref_peak
diff_final = 100 * abs(final_power - ref_final) / ref_final
plt.semilogy(tsteps, diff_peak, 'kx-')
plt.semilogy(tsteps, diff_final, 'rx-')
plt.xlabel('Number of Time Steps')
plt.ylabel('Error in Power (%)')
plt.legend(['Peak Power', 'Final Power'])
plt.show()

# find number of steps to converge to 1% peak power
peak_steps = tsteps[-1]
for i, error in enumerate(diff_peak):
    if error < 10:
        peak_steps = tsteps[i]
        break;
print "peak power convergence delta-t = ", 2.0/(peak_steps-1)

# find number of steps to converge to 1% final power
for i, error in enumerate(diff_final):
    if error < 1:
        final_steps = tsteps[i]
        break;
print "final power convergence delta-t = ", 2.0/(final_steps-1)

# calculate convergence for problem requested time steps
powers = range(7)
mult = 2**np.array(powers)
peak_power = []
final_power = []
for i, m in enumerate(mult):

    npts = (peak_steps-1)/m + 1
    print "Solving transient with ", npts, " time steps"
    
    # set time steps
    t = np.linspace(0,2,npts)
    trans.setCalcTimes(t)

    # solve reactor kinetics problem
    result = rk.solveTransient(trans, rkParams)
    P = np.array(result.getPower())
    P = P/P[0]

    # record power peak
    peak_power.append(max(P))
    final_power.append(P[-1])

    # plot converged core power vs time
    if i == 0:
        plt.plot(t, P/P[0], 'k.-')
        plt.xlabel('Time (s)')
        plt.ylabel('Relative Power')
        plt.show()

#FIXME
ref_peak = peak_power[0]
ref_final = final_power[0]

# plot error in peak power and final power
tsteps = 2.0 / ( (peak_steps-1) / np.array(mult) )
peak_power = np.array(peak_power)
final_power = np.array(final_power)
diff_peak = 100 * abs(peak_power - ref_peak) / ref_peak
diff_final = 100 * abs(final_power - ref_final) / ref_final
plt.semilogy(tsteps, diff_peak, 'kx-')
plt.semilogy(tsteps, diff_final, 'rx-')
plt.xlabel('Time Step Size (s)')
plt.ylabel('Error in Power (%)')
plt.legend(['Peak Power', 'Final Power'], loc=5)
plt.show()

# evaluate fractional error in peak, final powers vs spatial mesh
spacing = [1,2,5,10]

peak_power = []
final_power = []
for dx in spacing:

    print "Solving transient for ", dx, " cm mesh spacing"
    
    # form all new mesh
    npts = [node/dx for node in nodes]
    mesh = []
    for i in range(2):
        temp = rk.Mesh()
        temp.setMesh(npts)
        temp.setWidths(widths)
        temp.setMaterials(mats[i])
        temp.setBoundaryConditions(2,2)
        mesh.append(temp)

    # set transient
    trans = rk.Transient()
    trans_t = [0, 0.1, 0.2, 0.5, 1.5, 2.0]
    trans_m = [1, 1, 0, 0, 1, 1]
    trans.setInterpTimes(trans_t)
    meshVector = []
    for i in xrange(len(trans_t)):
        meshVector.append(mesh[ trans_m[i] ])
    trans.setMeshVector(meshVector)

    # set time steps
    t = np.linspace(0, 2, peak_steps)
    trans.setCalcTimes(t)

    # solve reactor kinetics problem
    result = rk.solveTransient(trans, rkParams)
    P = np.array(result.getPower())
    P = P/P[0]

    # record power peak, final power
    peak_power.append(max(P))
    final_power.append(P[-1])

# plot errors in peak, final power
ref_peak = peak_power[0]
ref_final = final_power[0]
peak_power = np.array(peak_power)
final_power = np.array(final_power)
diff_peak = 100 * abs(peak_power - ref_peak) / ref_peak
diff_final = 100 * abs(final_power - ref_final) / ref_final
plt.plot(spacing, diff_peak, 'kx-')
plt.plot(spacing, diff_final, 'rx-')
plt.xlabel('Mesh Spacing (cm)')
plt.ylabel('Error in Power (%)')
plt.legend(['Peak Power', 'Final Power'])
plt.show()
'''
