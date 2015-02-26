#include"diffSolver.h"

/*
    Forms the loss and fission matricies from given XSdata and a defined mesh
    with boundary conditions BC
    Arguments:
        mesh:   Mesh struct that contains region spacing and cross sections
        BC:     BoundaryConditions struct that contains values describing the
                desired boundary conditions associated with the problem
                options:
                0 = zero flux BC
                1 = zero incoming flux BC
                2 = reflective BC
        A:      Loss matrix
        F:      Fission matrix
        index:  Indexing function that calculates the indecies into matricies
                and vectors.
*/
void formMatrixProblem(Mesh mesh, BoundaryConditions BC, 
        Sparse &A, Sparse &F, Indexer index)
{
   
    // extract dimensions 
    int N = index._N;
    int G = index._G;

    // form Dhat values
    std::vector<std::vector<double> > Dhat;
    for(int n=0; n<N-1; n++)
    {
        /* load material and mesh spacing of current and 
           neighboring cells */
        double delta1 = mesh.delta[n];
        double delta2 = mesh.delta[n+1];
        XSdata* mat1 = mesh.material[n];
        XSdata* mat2 = mesh.material[n+1];

        // calculate Dhat for all groups
        std::vector<double> temp;
        for(int g=0; g<G; g++)
        {
            // load diffusion coefficients for group g
            double D1 = mat1->D[g];
            double D2 = mat2->D[g];
            
            // calculate Dhat and add to vector
            double dhat_value = 2*D1*D2 / (delta2*D1 + delta1*D2);
            temp.push_back(dhat_value);
        }
        Dhat.push_back(temp);
    }
        
    // setup matricies
    double value;
    for(int n=0; n<N; n++)
    {
        double delta = mesh.delta[n];
        XSdata* mat = mesh.material[n];
        for(int g=0; g<G; g++)
        {
            // calculate removal term
            double diag_val;
            double sigr = mat->siga[g];
            for(int gp=0; gp<G; gp++)
            {
                if(gp != g)
                    sigr += mat->sigs[g][gp];
            }

            // calculate loss matrix terms
            if(n == 0 or n == N-1)
            {

                // treat boundary conditions
                double Dhat_neighbor;
                int bc;

                // add neighbor term
                if(n == 0)
                {
                    Dhat_neighbor = Dhat[n][g];
                    if(Dhat_neighbor != 0)
                        A.setVal(index(n,g), index(n+1,g), -Dhat_neighbor);
                    bc = BC.left;
                }
                else
                {
                    Dhat_neighbor = Dhat[n-1][g];
                    if(Dhat_neighbor != 0)
                        A.setVal(index(n,g), index(n-1,g), -Dhat_neighbor);
                    bc = BC.right;
                }

                // add diagonal term
                if(bc == 0)
                {
                    // zero flux BC
                    diag_val = 
                        Dhat_neighbor + 2*mat->D[g] / delta + sigr*delta;
                } 
                else if(bc ==1)
                {   
                    // zero incoming flux BC
                    diag_val = Dhat_neighbor + 
                        2*mat->D[g] / (delta + 4*mat->D[g]) + sigr*delta;
                }
                else if(bc == 2)
                {
                    // reflected BC
                    diag_val = Dhat_neighbor + sigr * delta;
                }
                else
                    std::printf("Error: unexpected boundary condition\n");

                A.setVal(index(n,g), index(n,g), diag_val);
            }
            else
            {
                // treat bulk terms

                // add neighbor terms
                A.setVal(index(n,g), index(n-1,g), -Dhat[n-1][g]);
                A.setVal(index(n,g), index(n+1,g), -Dhat[n][g]);

                // add diagonal elements
                diag_val = Dhat[n-1][g] + Dhat[n][g] + sigr*delta;
                A.setVal(index(n,g), index(n,g), diag_val);
            }

            // add scattering source terms
            for(int gp=0; gp<G; gp++)
            {
                if(gp != g)
                {
                    double sigs = mat->sigs[gp][g];
                    if(sigs != 0)
                        A.setVal(index(n,g), index(n,gp), -sigs*delta);
                }
            }

            // calculate fission matrix terms
            for(int gp=0; gp<G; gp++)
            {
                double fission = mat->chi[g]*mat->nuSigf[gp]*delta;
                if(fission != 0)
                    F.setVal(index(n,g), index(n,gp), fission);
            }
        }
    }

    return;
}

/*
Fills geometry with materials and forms a mesh
Arguments:
    materials:  list of materials to fill each region
    widths:     list of region widths
    npts:       number of descritization points in each region

Output:
    mesh:   a Mesh object that contains XS information and cell widths  
*/
    
Mesh fillGeometry(std::vector<XSdata*> materials, std::vector<double> widths, 
        std::vector<int> npts)
{    
    // extract the total number of regions
    int num_regions = materials.size();

    // calculate total mesh points
    int N = 0;
    for(int n=0; n<npts.size(); n++)
        N += npts[n];

    // initialize mesh
    Mesh mesh;
    mesh.delta = std::vector<double> (N,0);
    
    // determine the number of groups
    int G = materials[0]->D.size();

    // form mesh and fill with materials
    int c = 0;
    for(int i=0; i<num_regions; i++)
    {
        for(int j=0; j<npts[i]; j++)
        {
            // fill XS data
            mesh.material.push_back( materials[i] );
            
            // calculate mesh
            mesh.delta[c] = (double) widths[i] / (double) npts[i];

            // increment counter
            c += 1;
        }
    }
    return mesh;
}

/*
    Solves the eigenvalue problem Ax = 1/k Fx <-> inv(A)Fx = kx
    Arguments:
        - A:            loss matrix
        - F:            fission matrix
        - N:            number of nodes
        - G:            number of energy groups
        - tol:          fission source iteration tolerance
        - inner_tol:    flux solver tolerance
        - maxiters:     maximum fission source iterations 
        - inner_solver: an integer identifying the flux solution solver
                        options: 0 = Point Jacobi, 1 = Gauss-Seidel

    Output:
        An dictionary containing solution information with keys:
            - 'flux': solution vector of the eigenvalue problem
            - 'keff': eigenvalue of the solution (largest eigenvalue)
            - 'power': the fission production corresponding to the solution
            - 'outer_iters': number of fission source iterations needed
            - 'inner_iters': average number of flux solution iterations needed
*/
eigenSolution eigen_solver(
        Sparse A, 
        Sparse F, 
        int N, 
        int G, 
        double tol, 
        double inner_tol, 
        int maxiters, 
        int inner_solver)
{
    // intialize flux to be flat, critical
    std::vector<double> phi ( A.size(), 1.0/double(N) );
    std::vector<double> phi_new ( A.size(), 0 );
    double k = 1;

    // allocate source
    std::vector<double> S;

    // track number of inner iterations
    int sum_inner_iters = 0;
    int outer_iters = 0;

    // outer iterations to converge eigenvalue problem
    for(int i=0; i < maxiters; i++)
    {
        // compute source
        S = F * phi;
        
        // solve inner iterations
        if(inner_solver == 0)
            phi_new = pointJacobi(A, S, inner_tol, maxiters, sum_inner_iters);
        
        else if(inner_solver == 1)
            phi_new = gaussSeidel(A, S, inner_tol, maxiters, sum_inner_iters);
        
        else
            std::cout << "Error: inner solver not recognized" << endl;

        // compute new eigenvalue (rayleigh quotient)
        k = dot(phi_new, phi) / dot(phi, phi);

        // create new fission sources
        std::vector<double> temp = F * phi_new;
        
        // calculate residual
        double res = 0;
        for(int n=0; n < temp.size(); n++)
        {
            if(temp[n] != 0)
                res += pow(k*S[n] / temp[n] - 1, 2);
        }
        res /= N;
        res = sqrt(res);

        // calculate new total fission source
        double sum_fission = 0;
        for(int n=0; n<temp.size(); n++)
            sum_fission += temp[n];

        // normalize phi
        for(int n=0; n<phi.size(); n++)
            phi[n] = phi_new[n] * N * G / sum_fission;


        // check for convergence
        std::cout << "res = " << res << ", iter# = " << sum_inner_iters/(i+1) << endl;
        if(res < tol)
        {
            outer_iters = i+1;
            break;
        }
    }

    // check if the maximum iterations exceeded
    if(outer_iters == 0)
    {
        printf("Warning: Maximum fission source iterations exceeded!\n");
        outer_iters = maxiters;
    }

    // structure eigenvalue result
    eigenSolution result;
    result.flux = phi;
    result.keff = k;
    result.outer_iters = outer_iters;
    result.power = F*phi;
    result.inner_iters = sum_inner_iters / (outer_iters);
    
    std::cout << "Source converged in " << outer_iters << 
        " outer iterations" << endl;
    return result;
}

/*
    Gauss-Seidel method to solve Ax = b
    Arguments:
        - A:        matrix A
        - b:        vector b
        - tol:      convergence tolerance
        - maxiters: maximum iterations
        - sumiters: number of iterations needed to converge problem
    Returns:
        The solution vector x
    */
std::vector<double> gaussSeidel(
        Sparse A, 
        std::vector<double> b, 
        double tol, 
        int maxiters, 
        int &sum_iters)
{
    // initialize solution
    int N = b.size();
    std::vector<double> x (N, 1);
    int iters = 0;

    // cycle through iterations
    for(int n=0; n<maxiters; n++)
    {
        // create new vector xnew which is a copy of x
        std::vector<double> xnew (x.size(), 0);

        // iterate through rows of solution vector
        for(int i=0; i < N; i++)
        {
            // initialize to rhs value
            xnew[i] = b[i];

            // iterate through columns of matrix 
            for( int j=0; j < N; j++)
            {
                // subtract upper traingular
                if( i < j )
                    xnew[i] -= A(i,j) * x[j];
                
                // subtract lower triangular (updated)
                else if( i > j )
                    xnew[i] -= A(i,j) * xnew[j];
            }

            // divide by diagonal
            xnew[i] /= A(i,i);
        }

        // copmpute difference between previous iteration
        std::vector<double> diff (x.size(), 0);
        for(int i=0; i < diff.size(); i++)
            diff[i] = (xnew[i] - x[i]) / x[i];
        double res = sqrt( dot(diff,diff) / N);
        
        // update solution vector x
        x = xnew;

        // check convergence
        if(res < tol)
        {
            iters = n+1;
            break;
        }
    }
    
    // check of maxiters exceeded
    if( iters == 0 )
    {
        std::cout << "Warning: Maximum inner iterations exceeded!" << endl;
        iters = maxiters;
    }

    // add to iteration count
    sum_iters += iters;
    
    return x;
}

/*
    Point Jacobi method to solve Ax = b
    Arguments:
        - A:        matrix A
        - b:        vector b
        - tol:      convergence tolerance
        - maxiters: maximum iterations
        - sumiters: number of iterations needed to converge problem
    Returns:
        The solution vector x
    */
std::vector<double> pointJacobi(
        Sparse A, 
        std::vector<double> b, 
        double tol, 
        int maxiters, 
        int &sum_iters)
{
    // initialize solution
    int N = b.size();
    std::vector<double> x (N, 1);
    int iters = 0;

    // cycle through iterations
    for(int n=0; n<maxiters; n++)
    {
        // create new vector xnew which is a copy of x
        std::vector<double> xnew (x.size(), 0);

        // iterate through rows of solution vector
        for(int i=0; i < N; i++)
        {
            // initialize to rhs value
            xnew[i] = b[i];

            // iterate through columns of matrix 
            for( int j=0; j < N; j++)
            {
                // subtract upper traingular
                if( i < j )
                    xnew[i] -= A(i,j) * x[j];
                
                // subtract lower triangular
                else if( i > j )
                    xnew[i] -= A(i,j) * x[j];
            }

            // divide by diagonal
            xnew[i] /= A(i,i);
        }

        // copmpute difference between previous iteration
        std::vector<double> diff (x.size(), 0);
        for(int i=0; i < diff.size(); i++)
            diff[i] = (xnew[i] - x[i]) / x[i];
        double res = sqrt( dot(diff,diff) / N);
        
        // update solution vector x
        x = xnew;

        // check convergence
        if(res < tol)
        {
            iters = n+1;
            break;
        }
    }
    
    // check of maxiters exceeded
    if( iters == 0 )
    {
        std::cout << "Warning: Maximum inner iterations exceeded!" << endl;
        iters = maxiters;
    }

    // add to iteration count
    sum_iters += iters;
    
    return x;
}



