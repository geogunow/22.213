from diff_solver import *
from utils import *

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size' : 18}
plt.rc('font', **font)

'''
script for pset 2
'''
# Define cross sections / nuclear data

# define simple homogeneous XS for testing
homo_mat = dict()
homo_mat['D'] = [1, 1]
homo_mat['siga'] = [1, 1]
homo_mat['sigs'] = [[1, 1],[0, 1]]
homo_mat['nuSigf'] = [1, 1]
homo_mat['chi'] = [1, 0]

# define water material
water = dict()
water['D'] = [1.5, 0.5]
water['siga'] = [0.0002, 0.010]
water['sigs'] = [[0, 0.032],[0, 0]]
water['nuSigf'] = [0, 0]
water['chi'] = [1, 0]

# define fuel1 material (3% enriched)
fuel1 = dict()
fuel1['D'] = [1.3, 0.5]
fuel1['siga'] = [0.0098, 0.114]
fuel1['sigs'] = [[0, 0.022],[0, 0]]
fuel1['nuSigf'] = [0.006, 0.195]
fuel1['chi'] = [1, 0]

# define fuel2 material (4% enriched)
fuel2 = dict()
fuel2['D'] = [1.3, 0.5]
fuel2['siga'] = [0.0098, 0.114]
fuel2['sigs'] = [[0, 0.022],[0, 0]]
fuel2['nuSigf'] = [0.008, 0.23]
fuel2['chi'] = [1, 0]

# define black absorber material
black = dict()
black['D'] = [1.3, 0.5]
black['siga'] = [999.99, 999.9]
black['sigs'] = [[0, 0.020],[0, 0]]
black['nuSigf'] = [0, 0]
black['chi'] = [1, 0]

# define fuel3 material (3% enriched + rod)
fuel3 = dict()
fuel3['D'] = [1.3, 0.5]
fuel3['siga'] = [0.0098, 0.118]
fuel3['sigs'] = [[0, 0.022],[0, 0]]
fuel3['nuSigf'] = [0.006, 0.195]
fuel3['chi'] = [1, 0]


# define physical dimensions & boudary conditions
BC = [2, 2]
materials = [water, fuel1, fuel3, fuel1, water]
widths = [30.0, 170.0, 20.0, 170.0, 30.0]
nodes = [3, 17, 2, 17, 3]
N = sum(nodes)

'''
Parts B & C (and implicitly A)
'''
# define mesh refinement
refinements = range(1,30)

# create arrays for tracking l2 norms and 
nodal_power = []
keff = np.zeros( len(refinements) )
dominance_ratio = np.zeros( len(refinements) )

# cycle through mesh refinements
for i, refine in enumerate(refinements):
    
    print "Mesh refinement ", refine
    
    # define the mesh
    npts = [pt*refine for pt in nodes]

    # form mesh and XS data
    XSdata, delta = fillGeometry(materials, widths, npts)

    # form loss and fission matricies
    A, F = formMatrixProblem(XSdata, delta, BC)

    # solve eigenvalue problem
    result = eigen_solver(A, F, k=1, N = N)
    keff[i] = result['keff']
    nodal_avg_power = calculate_nodal_avg(result['power'], refine)
    nodal_power.append( nodal_avg_power )

    # compute dominance ratio
    dominance_ratio[i] = compute_dominance_ratio(A, F)

    # plot fluxes
    if refine == refinements[-1]:
        plot_fluxes(delta, result['flux'], sum(npts))

    print "Keff = ", keff[i]

# calculate errors relative to finest mesh
keff_errors = calculate_scalar_errors(keff[:-1], keff[-1])*10**5
nodal_power_errors = calculate_l2_errors(nodal_power[:-1], nodal_power[-1], N)

# plot error trends
dual_plot(refinements[:-1], keff_errors, nodal_power_errors, 
        xlabel='Mesh Points per Node', ylabel1='Relative k-eff Error (pcm)',
        ylabel2='Nodal Power RMS Error')

# plot error trends
dual_log_plot(refinements[:-1], keff_errors, nodal_power_errors, 
        xlabel='Mesh Points per Node', ylabel1='Relative k-eff Error (pcm)',
        ylabel2='Nodal Power RMS Error')

# plot dominance ratio
plt.semilogy(refinements, dominance_ratio, 'kx-')
plt.xlabel('Mesh Points per Node')
plt.ylabel('Dominance Ratio')
plt.show()

# plot dominance ratio
plt.plot(refinements, dominance_ratio, 'kx-')
plt.xlabel('Mesh Points per Node')
plt.ylabel('Dominance Ratio')
plt.show()

'''
Parts D and E
'''
# form mesh and XS data
XSdata, delta = fillGeometry(materials, widths, nodes)

# form loss and fission matricies
A, F = formMatrixProblem(XSdata, delta, BC)

# define toleratnces for inner solver
powers = range(-5,0)
tolerances = [10**p for p in powers]

# determine behavior of point jacobi and gauss-seidel
for solver in ['PJ', 'GS']:

    # initialize vector of iteration numbers
    outer_iters = []
    inner_iters = []

    for inner_tol in tolerances:
        
        # solve eigenvalue problem
        result = eigen_solver(A, F, k=1, N = N, inner_solver=solver, 
                inner_tol=inner_tol)

        # save iteration results
        outer_iters.append( result['outer_iters'] )
        inner_iters.append( result['inner_iters'] )

    # plot iteration results
    xlabel = 'Tolerance on L-2 norm of Flux Solutions' 
    ylabel1 = 'Flux Solution Iteration Count'    
    ylabel2 = 'Fission Source Iteration Count'

    dual_logx_plot(tolerances, inner_iters, outer_iters, xlabel=xlabel,
            ylabel1=ylabel1, ylabel2=ylabel2)

'''
Part F
'''
# refine mesh
npts = [pt*20 for pt in nodes]

# form mesh and XS data
XSdata, delta = fillGeometry(materials, widths, npts)

# form loss and fission matricies
A, F = formMatrixProblem(XSdata, delta, BC)
A = A.getT()
F = F.getT()
result = eigen_solver(A, F, k=1, N = N)

# plot fluxes and print k
plot_fluxes(delta, result['flux'], sum(npts))
print "k = ", result['keff']
