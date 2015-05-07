#include"diffSolver.h"

/*
   Function for solving a diffusion problem given a mesh with boundary
   conditions. This function is designed for easy communication with Python
   via Swig.
    Arguments:
        - mesh:         Mesh object that contains region spacing, cross 
                        sections, and boundary condition information
        - outer_tol:    fission source iteration tolerance
        - inner_tol:    flux solver tolerance
        - maxiters:     maximum fission source iterations 
        - inner_solver: an integer identifying the flux solution solver
                        options: 0 = Point Jacobi, 1 = Gauss-Seidel

    Output:
        An eigenSolution struct containing the following information:
            - 'flux': solution vector of the eigenvalue problem
            - 'keff': eigenvalue of the solution (largest eigenvalue)
            - 'power': the fission production corresponding to the solution
            - 'outer_iters': number of fission source iterations needed
            - 'inner_iters': average number of flux solution iterations needed
 */
eigenSolution solveCritical(Mesh mesh, double outer_tol, double inner_tol,
        int maxiters, int inner_solver)
{
    // error checking
    eigenSolution NULL_solution = eigenSolution();
    if( mesh._N == 0 )
    {
        std::cout << "Error: failed to set mesh correctly" << endl;
        return NULL_solution;
    }
    if( mesh._G == 0)
    {
        std::cout << "Error: failed to set materials correctly" << endl;
        return NULL_solution;
    }
    if( mesh._BC.left == -1)
    {
        std::cout << "Error: failed to set boundary conditions" << endl;
        return NULL_solution;
    }
    if( mesh._BC.right == -1)
    {
        std::cout << "Error: failed to set boundary conditions" << endl;
        return NULL_solution;
    }
    if( mesh._delta.size() == 0)
    {
        std::cout << "Error: failed to set mesh widths" << endl;
        return NULL_solution;
    }

    // initialize matrices
    int N = mesh._N;
    int G = mesh._G; 
    Indexer index = Indexer(N,G);
   
    // form matrices
    Sparse A = Sparse(N*G, N*G);
    Sparse F = Sparse(N*G, N*G);
    formSteadyStateMatrixProblem(mesh, A, F, index);
        
    // solve eigenvalue problem
    N = mesh._nodes;
    eigenSolution result = eigen_solver(A, F, N, G, outer_tol, inner_tol, 
            maxiters, inner_solver);

    // index soltion vectors by groups
    result.indexArrays(index);

    return result;
}

// form SigA matrix
Sparse formSigAMatrix(Mesh mesh, Indexer index)
{
    // extract dimensions 
    int N = index._N;
    int G = index._G;

    // intialize sparse matrix
    Sparse SigA = Sparse(N*G, N*G);
    
    // setup matricies
    for(int n=0; n<N; n++)
    {
        double delta = mesh._delta[n];
        XSdata* mat = mesh._material[n];
        for(int g=0; g<G; g++)
        {
            // calculate removal term
            SigA.setVal( index(n,g), index(n,g), mat->siga[g]*delta);
        }
    }
    return SigA;
}

// form SigS matrix
Sparse formSigSMatrix(Mesh mesh, Indexer index)
{
   
    // extract dimensions 
    int N = index._N;
    int G = index._G;

    // initialize matrix
    Sparse SigS = Sparse(N*G, N*G);

       
    // setup matrix
    for(int n=0; n<N; n++)
    {
        double delta = mesh._delta[n];
        XSdata* mat = mesh._material[n];
        for(int g=0; g<G; g++)
        {
            // calculate removal term (without absorption)
            double sigr = 0;
            for(int gp=0; gp<G; gp++)
            {
                if(gp != g)
                    sigr += mat->sigs[g][gp];
            }
            SigS.setVal(index(n,g), index(n,g), sigr*delta);
                        
            // add scattering source terms
            for(int gp=0; gp<G; gp++)
            {
                if(gp != g)
                {
                    double sigs = mat->sigs[gp][g];
                    if(sigs != 0)
                        SigS.setVal(index(n,g), index(n,gp), -sigs*delta);
                }
            }
        }
    }
    return SigS;
}

// form Dhat matrix
Sparse formDhatMatrix(Mesh mesh, Indexer index)
{
    // extract dimensions
    int N = index._N;
    int G = index._G;

    // intialize matrix
    Sparse Dhat = Sparse(N*G, N*G);

    // form Dhat values
    std::vector<std::vector<double> > Dhat_vals;
    for(int n=0; n<N-1; n++)
    {
        /* load material and mesh spacing of current and 
           neighboring cells */
        double delta1 = mesh._delta[n];
        double delta2 = mesh._delta[n+1];
        XSdata* mat1 = mesh._material[n];
        XSdata* mat2 = mesh._material[n+1];

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
        Dhat_vals.push_back(temp);
    }
        
    // setup matrix
    for(int n=0; n<N; n++)
    {
        double delta = mesh._delta[n];
        XSdata* mat = mesh._material[n];
        for(int g=0; g<G; g++)
        {
            // calculate removal term
            double diag_val;
            
            // calculate loss matrix terms
            if(n == 0 || n == N-1)
            {
                // treat boundary conditions
                double Dhat_neighbor;
                int bc;

                // add neighbor term
                if(n == 0)
                {
                    Dhat_neighbor = Dhat_vals[n][g];
                    if(Dhat_neighbor != 0)
                        Dhat.setVal(index(n,g), index(n+1,g), -Dhat_neighbor);
                    bc = mesh._BC.left;
                }
                else
                {
                    Dhat_neighbor = Dhat_vals[n-1][g];
                    if(Dhat_neighbor != 0)
                        Dhat.setVal(index(n,g), index(n-1,g), -Dhat_neighbor);
                    bc = mesh._BC.right;
                }

                // add diagonal term
                if(bc == 0)
                {
                    // zero flux BC
                    diag_val = Dhat_neighbor + 2*mat->D[g] / delta;
                } 
                else if(bc ==1)
                {   
                    // zero incoming flux BC
                    diag_val = Dhat_neighbor + 
                        2*mat->D[g] / (delta + 4*mat->D[g]);
                }
                else if(bc == 2)
                {
                    // reflected BC
                    diag_val = Dhat_neighbor;
                }
                else
                    std::printf("Error: unexpected boundary condition\n");

                Dhat.setVal(index(n,g), index(n,g), diag_val);
            }
            else
            {
                // treat bulk terms

                // add neighbor terms
                Dhat.setVal(index(n,g), index(n-1,g), -Dhat_vals[n-1][g]);
                Dhat.setVal(index(n,g), index(n+1,g), -Dhat_vals[n][g]);

                // add diagonal elements
                diag_val = Dhat_vals[n-1][g] + Dhat_vals[n][g];
                Dhat.setVal(index(n,g), index(n,g), diag_val);
            }
        }
    }
    return Dhat;
}


// form nuSigf matrix
Sparse formFMatrix(Mesh mesh, Indexer index)
{
    // extract dimensions
    int N = index._N;
    int G = index._G;

    // intialize matrix
    Sparse F = Sparse(N*G, N*G);
    
    // setup matrix
    for(int n=0; n<N; n++)
    {
        double delta = mesh._delta[n];
        XSdata* mat = mesh._material[n];
        for(int g=0; g<G; g++)
        {
            // calculate fission matrix terms
            for(int gp=0; gp<G; gp++)
            {
                double fission = mat->chi[g]*mat->nuSigf[gp]*delta;
                if(fission != 0)
                    F.setVal(index(n,g), index(n,gp), fission);
            }
        }
    }

    return F;
}   



/*
    Forms the loss and fission matricies from given XSdata and a defined mesh
    with boundary conditions BC
    Arguments:
        mesh:   Mesh object that contains region spacing and cross sections
                0 = zero flux BC
                1 = zero incoming flux BC
                2 = reflective BC
        A:      Loss matrix
        F:      Fission matrix
        index:  Indexing function that calculates the indecies into matricies
                and vectors.
*/
void formSteadyStateMatrixProblem(Mesh mesh, Sparse &A, Sparse &F, 
        Indexer index)
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
        double delta1 = mesh._delta[n];
        double delta2 = mesh._delta[n+1];
        XSdata* mat1 = mesh._material[n];
        XSdata* mat2 = mesh._material[n+1];

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
    for(int n=0; n<N; n++)
    {
        double delta = mesh._delta[n];
        XSdata* mat = mesh._material[n];
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
            if(n == 0 || n == N-1)
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
                    bc = mesh._BC.left;
                }
                else
                {
                    Dhat_neighbor = Dhat[n-1][g];
                    if(Dhat_neighbor != 0)
                        A.setVal(index(n,g), index(n-1,g), -Dhat_neighbor);
                    bc = mesh._BC.right;
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
        An eigenSolution struct containing the following information:
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
            phi_new = A.pointJacobi(S, phi, inner_tol, maxiters, sum_inner_iters);
        
        else if(inner_solver == 1)
            phi_new = A.gaussSeidel(S, phi, inner_tol, maxiters, sum_inner_iters);
        
        else if(inner_solver == 2)
            phi_new = A.optimalSOR(S, phi, inner_tol, maxiters, sum_inner_iters);
        
        else
            std::cout << "Error: inner solver not recognized" << endl;

        // compute new eigenvalue (rayleigh quotient)
        k = dot(phi_new, phi) / dot(phi, phi);

        // create new fission sources
        std::vector<double> temp = F * phi_new;
        
        // calculate new total fission source
        double sum_fission = 0;
        for(int n=0; n<temp.size(); n++)
            sum_fission += temp[n];
        
        // calculate residual
        double res = 0;
        for(int n=0; n < temp.size(); n++)
        {
            if(temp[n] != 0)
                res += pow(k*S[n] / temp[n] - 1, 2);
        }
        res /= (N);
        res = sqrt(res);

        // normalize phi
        for(int n=0; n<phi.size(); n++)
            phi[n] = phi_new[n] * N * G / sum_fission;

        // check for convergence
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
    eigenSolution result = eigenSolution();
    result.flux = phi;
    result.keff = k;
    result.outer_iters = outer_iters;
    result.power = F*phi;
    result.inner_iters = sum_inner_iters / (outer_iters);
    
    return result;
}
