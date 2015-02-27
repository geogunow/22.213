import rksolver as rk
import numpy as np
import matplotlib.pyplot as plt
from utils import *

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

# set cross sections for fuel3 (3% enriched + rod)
fuel3 = rk.XSdata()
fuel3.setD( [1.3, 0.5] )
fuel3.setSiga( [0.0098, 0.118] )
fuel3.setSigs( [[0, 0.022], [0, 0]])
fuel3.setNuSigf( [0.006, 0.195] )
fuel3.setChi( [1.0, 0.0] )

# describe geometry
mats = [water, fuel1, fuel3, fuel1, water]
widths = [30.0, 170.0, 20.0, 170.0, 30.0]
nodes = [3, 17, 2, 17, 3]
N = sum(nodes)

# define mesh refinement
refinements = range(1,12)

# create arrays for tracking l2 norms and keff
nodal_power = []
keff = np.zeros( len(refinements) )

# eigenvalue solver parameters
outer_tol = 10**-6
inner_tol = 10**-5
maxiters = 10000
inner_solver = 1 # optimal SOR

'''
Parts A and B
'''
# cycle through mesh refinements
for i, refine in enumerate(refinements):
    
    print "Mesh refinement ", refine

    # define the mesh
    npts = [pt*refine for pt in nodes]
    mesh = rk.Mesh()
    mesh.setMesh(npts)
    mesh.setWidths(widths)
    mesh.setMaterials(mats)
    mesh.setBoundaryConditions(2,2)
    x = mesh.getX()

    # solve the eigenvalue problem
    solution = rk.solveDiffusion(mesh, outer_tol, inner_tol, maxiters,
            inner_solver)

    # save power
    total_power_array = []
    for g in range(G):
        groupPower = solution.getPower(g+1)
        nodal_avg_power = calculate_nodal_avg(groupPower, refine)
        total_power_array += list(nodal_avg_power)
    nodal_power.append(total_power_array)

    # save keff
    keff[i] = solution.keff

    # plot flux
    if refine == refinements[-1]:
        for g in range(G):
            plt.plot(x, solution.getFlux(g+1), 'x-')
        plt.xlabel("x (cm)")
        plt.ylabel("Normalized Flux")
        plt.show()

    print "Keff = ", solution.keff

# calculate RMS errors relative to finest mesh
keff_errors = calculate_scalar_errors(keff[:-1], keff[-1])*10**5
nodal_power_errors = calculate_RMS_errors(nodal_power[:-1], nodal_power[-1], N)

# plot error trends
dual_plot(refinements[:-1], keff_errors, nodal_power_errors,
        xlabel='Mesh Points per Node', ylabel1='Relative k-eff Error (pcm)',
        ylabel2='Nodal Power RMS Error')

# plot error trends
dual_log_plot(refinements[:-1], keff_errors, nodal_power_errors,
        xlabel='Mesh Points per Node', ylabel1='Relative k-eff Error (pcm)',
        ylabel2='Nodal Power RMS Error')

'''
Parts D and E
'''
# form mesh
mesh = rk.Mesh()
mesh.setMesh(nodes)
mesh.setWidths(widths)
mesh.setMaterials(mats)
mesh.setBoundaryConditions(2,2)

# define tolerances for inner solver
powers = range(-5,0)
tolerances = [10**p for p in powers]
 
# determine behavior of point jacobi and gauss-seidel
for solver in [0,1]:

    # set linear solver
    inner_solver = solver

    # initialize vector of iteration numbers
    outer_iters = []
    inner_iters = []

    for inner_tol in tolerances:

        # solve eigenvalue problem
        solution = rk.solveDiffusion(mesh, outer_tol, inner_tol, maxiters,
                inner_solver)

        # save iteration results
        outer_iters.append(solution.outer_iters)
        inner_iters.append(solution.inner_iters)

    # plot iteration results
    xlabel = 'Tolerance on RMS norm of Flux Solutions'
    ylabel1 = 'Flux SOlution Iteration Count'
    ylabel2 = 'Fission Source Iteration Count'

    dual_logx_plot(tolerances, inner_iters, outer_iters, xlabel=xlabel,
            ylabel1=ylabel1, ylabel2=ylabel2)

'''
Part F
'''
#x = solution.getFlux(1)
