#include"diff_solver.h"

using namespace std;

/*
Forms the loss and fission matricies from given XSdata and a defined mesh with
boundary conditions BC
Arguments:
    XSdata: cross section data for each mesh point formed by the fillGeometry
            function
    delta:  the mesh of the problem corresponding to XSdata
    BC:     an array of two values describing the desired boundary conditions 
            associated with the problem
            options:
            0 = zero flux BC
            1 = zero incoming flux BC
            2 = reflective BC
    A:      loss matrix
    F:      fission matrix
*/
void formMatrixProblem(map<string,double> XSdata, vector<double> delta, 
        vector<int> BC, smatrix A, smatrix F){
    
    // initialize matricies
    N = len(D)
    G = len(D[0])
    A = np.zeros( (N*G, N*G) )
    F = np.zeros( (N*G, N*G) )

    // define indexing function into matricies
    def index(n, g):
        #return G*n + g # groups inside
        return N*g + n # groups outside

    # form Dhat values
    Dhat = np.zeros( (N-1, G) )
    for n in xrange(N-1):
        for g in xrange(G):
            Dhat[n, g] = 2*D[n ,g]*D[n+1, g] \
                    / (delta[n+1] * D[n, g] + delta[n] * D[n+1, g])

    # setup matricies
    for g in xrange(G):
        for n in xrange(N):

            # calculate removal term
            sigr = siga[n, g]
            for gp in xrange(G):
                if gp != g:
                    sigr += sigs[n, g, gp]

            # calculate loss matrix terms
            if n == 0 or n == N-1:

                # treat boundary conditions

                # add neighbor term
                if n == 0:
                    Dhat_neighbor = Dhat[n,g]
                    A[index(n,g), index(n+1,g)] -= Dhat_neighbor
                    bc = BC[0]
                else:
                    Dhat_neighbor = Dhat[n-1,g]
                    A[index(n,g), index(n-1,g)] -= Dhat_neighbor
                    bc = BC[1]

                # add diagonal term
                if bc == 0:

                    # zero flux BC
                    A[index(n,g), index(n,g)] += Dhat_neighbor + \
                            2*D[n,g] / delta[n] + sigr * delta[n]
                
                elif bc == 1:

                    # zero incoming flux BC
                    A[index(n,g), index(n,g)] += Dhat_neighbor + \
                            2 * D[n,g] / (delta[n] + 4*D[n,g]) + sigr * delta[n]
                
                elif bc == 2:

                    # reflected BC
                    A[index(n,g), index(n,g)] += \
                            Dhat_neighbor + sigr * delta[n]

                else:
                    print "Error: unexpected boundary condition"
                    exit()

            else:

                # treat bulk terms

                # add neighbor terms
                A[index(n,g), index(n-1,g)] -= Dhat[n-1,g]
                A[index(n,g), index(n+1,g)] -= Dhat[n,g]

                # add diagonal elements
                A[index(n,g), index(n,g)] += Dhat[n-1,g] + Dhat[n,g] \
                        + sigr * delta[n]

            # add scattering source terms
            for gp in xrange(G):
                if gp != g:
                    A[index(n,g), index(n,gp)] -= sigs[n, gp, g] * delta[n]

            # calculate fission matrix terms
            for gp in xrange(G):
                F[index(n,g), index(n,gp)] += \
                        chi[n,g] * nuSigf[n,gp] * delta[n]

    return np.matrix(A), np.matrix(F)

'''
Fills geometry with materials and forms a mesh
Arguments:
    materials:  list of materials to fill each region
    widths:     list of region widths
    npts:       number of descritization points in each region

Output:
    XSdata: a dictionary containing the spatial distributions of 
            each cross section
    delta:  the mesh associated with the XSdata
'''
def fillGeometry(materials, widths, npts):
    
    # extract the total number of regions
    num_regions = len(materials)

    # calculate total mesh points
    N = 0
    for n in npts:
        N += n

    # initialize mesh
    delta = np.zeros(N)
    
    # determine the number of groups
    G = len( materials[0]['D'] )

    # initialize dictionaries of XS data
    XSdata = dict()
    attributes = ['D', 'siga', 'sigs', 'nuSigf', 'chi']
    for attr in attributes:
        if attr == 'sigs':
            XSdata[attr] = np.zeros( (N, G, G) )
        else:
            XSdata[attr] = np.zeros( (N, G) )

    # form mesh and fill with materials
    i = 0
    for (mat, w, n) in zip(materials, widths, npts):

        for j in xrange(n):
            
            # fill XS data
            for attr in attributes:
                XSdata[attr][i] = mat[attr]

            # calculate mesh
            delta[i] = float(w) / float(n)

            i += 1
        
    return XSdata, delta

'''
Solves the eigenvalue problem Ax = 1/k Fx <-> inv(A)Fx = kx
Arguments:
    - A:            loss matrix
    - F:            fission matrix
    - N:            number of nodes (optional)
                    default: None (recalculated to be number of mesh cells)
    - G:            number of energy groups (optional)
                    default: 2
    - tol:          fission source iteration tolerance (optional)
                    default: 1.0e-6
    - inner_tol:    flux solver tolerance (optional)
                    default: 1.0e-6
    - maxiters:     maximum fission source iterations (optional)
                    default: 1000
    - inner_solver: a string identifying the flux solution solver (optional)
                    default: intrinsic
                    options: GS, PJ, intrinsic

Output:
    An dictionary containing solution information with keys:
        - 'flux': solution vector of the eigenvalue problem
        - 'keff': eigenvalue of the solution (largest eigenvalue)
        - 'power': the fission production corresponding to the solution vector
        - 'outer_iters': number of fission source iterations needed
        - 'inner_iters': average number of flux solution iterations needed
'''
def eigen_solver(A, F, k=1, N=None, G=2, tol = 10**-6, inner_tol = 10**-6, 
        maxiters=1000, inner_solver='intrinsic'):

    # initialize solution
    ND = A.shape[0]
    if N == None:
        N = ND / G
    phi = np.ones( (ND,1) ) / N
    phi = np.matrix(phi)

    # initialize total fission source
    total_fission = F * phi

    # track number of inner iterations
    sum_inner_iters = 0

    # outer iterations to converge eigenvalue problem
    for i in xrange(maxiters):

        # compute source
        S = F * phi

        # solve inner iterations
        if inner_solver == 'GS':
            phi_new, inner_iters = gauss_seidel(A, S, tol=inner_tol)
            sum_inner_iters += inner_iters
        elif inner_solver == 'PJ':
            phi_new, inner_iters = point_jacobi(A, S, tol=inner_tol)
            sum_inner_iters += inner_iters
        else:
            phi_new = np.linalg.solve(A, S)

        # compute new eigenvalue (rayleigh quotient)
        k = np.dot(phi_new.getT(), phi) / np.dot(phi.getT(), phi)

        # create new fission sources
        total_fission = F * phi
        temp = F * phi_new
        
        # calculate residual
        res = 0
        for f_new, f_old in zip(temp, total_fission):
            if f_new != 0:
                res += (k*f_old / f_new - 1)**2
        res /= N
        res = sqrt(res)

        # normalize phi
        phi = phi_new * N * G / sum(temp)

        # check for convergence
        if res < tol:
            break

    # sturcture result and return answer
    result = dict()
    result['flux'] = phi
    result['keff'] = k
    result['outer_iters'] = i+1
    result['power'] = F*phi
    result['inner_iters'] = sum_inner_iters / (i+1)
    
    print "Source converged in ", result['outer_iters'], " outer iterations"
    return result

'''
Gauss-Seidel method to solve Ax = b
Arguments:
    - A:        matrix A
    - b:        vector b
    - maxiters: maximum iterations (optional)
                default: 1000
    - tol:      convergence tolerance
                default: 1.0e-6
Returns:
    A tuple containing:
        - the solution vector x
        - required iteration count
'''
def gauss_seidel(A, b, maxiters=1000, tol = 10**-6):
   
    N = len(b)
    x = np.ones( (N, 1) )
    x = np.matrix(x)

    # cycle through iterations
    for n in xrange(maxiters):

        # create new vector xnew
        xnew = copy(x)

        # iterate through rows of solution vector
        for i in xrange(N):

            # initialize with RHS value
            xnew[i] = b[i]

            # subtract upper traingular
            xnew[i] -= A[i,i+1:N] * x[i+1:N]

            # subtract lower triangular (updated)
            xnew[i] -= A[i,:i] * xnew[:i]

            # divide by diagonal
            xnew[i] /= A[i,i]

        # check for convergence
        diff = (xnew - x) / x
        res = sqrt(np.dot(diff.getT(), diff) / N)
        x = xnew
        if res < tol:
            break

    #print "Converged in ", n+1, " inner iterations"
    return x, n+1

'''
Point Jacobi method to solve Ax = b
Arguments:
    - A:        matrix A
    - b:        vector b
    - maxiters: maximum iterations (optional)
                default: 1000
    - tol:      convergence tolerance
                default: 1.0e-6
Returns:
    A tuple containing:
        - the solution vector x
        - required iteration count
'''
def point_jacobi(A, b, maxiters=1000, tol = 10**-6):
    
    N = len(b)
    x = np.ones( (N, 1) )
    x = np.matrix(x)

    # cycle through iterations
    for n in xrange(maxiters):

        # create new vector xnew
        xnew = copy(x)

        # iterate through rows of solution vector
        for i in xrange(N):

            # initialize with RHS value
            xnew[i] = b[i]

            # subtract upper traingular
            xnew[i] -= A[i,i+1:N] * x[i+1:N]

            # subtract lower triangular
            xnew[i] -= A[i,:i] * x[:i]

            # divide by diagonal
            xnew[i] /= A[i,i]

        # check for convergence
        diff = (xnew - x)/x
        res = sqrt(np.dot(diff.getT(), diff) / N)
        x = xnew
        if res < tol:
            break

    #print "Converged in ", n+1, " inner iterations"
    return x, n+1
