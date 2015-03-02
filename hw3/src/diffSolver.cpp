#include"diffSolver.h"

//TODO: write interface (may need to change to arrays)
void solveTransient(
        std::vector<Mesh> meshVector,
        std::vector<double> timeVector,
        std::vector<double> timeSteps,
        RKdata rkParams)
{
    // set transient variables
    double tol = pow(10,-6);
    int maxiters = pow(10,5);
    double outer_tol = pow(10,-6);
    double inner_tol = pow(10,-6);
    int inner_solver = 2;
    int sum_inner_iters = 0;

    //TODO: error checking
    std::vector<double> npower;


    // get first mesh to solve steady state
    Mesh mesh = meshVector[0];

    // initialize matrices
    Indexer index = Indexer(mesh._N,mesh._G);

    // set N to the number of nodes
    int N = mesh._nodes;
    
    // form steady state matrices
    Sparse F = formFMatrix(mesh, index);
    Sparse SigA = formSigAMatrix(mesh, index);
    Sparse SigS = formSigSMatrix(mesh, index);
    Sparse Dhat = formDhatMatrix(mesh, index);
    Sparse A = SigA + SigS + Dhat;
    
    // solve eigenvalue problem
    eigenSolution result = eigen_solver(A, F, N, mesh._G, outer_tol, inner_tol, 
            maxiters, inner_solver);
    
    // extract critical k for transient solution
    double kcrit = result.keff;

    // intialize precursor and flux structure
    std::vector<std::vector<std::vector<double> > > precursors;
    std::vector<std::vector<double> > flux;
    for(int t=0; t < timeSteps.size(); t++)
    {
        std::vector<std::vector<double> > precursor_1;
        std::vector<double> flux_1(mesh._N * mesh._G, 0);
        for(int n=0; n < mesh._N; n++)
        {
            std::vector<double> precursor_2 (rkParams.I, 0);
            precursor_1.push_back(precursor_2);
        }
        precursors.push_back(precursor_1);
        flux.push_back(flux_1);
    }

    // initialize first time step of flux
    flux[0] = result.flux;

    // calculate fission production in each cell
    std::vector<double> fissionProd(mesh._N, 0);
    for(int n=0; n < mesh._N; n++)
        for(int g=0; g < mesh._G; g++)
            fissionProd[n] += result.power[index(n,g)];
    
    // calculate intial precursor population
    for(int n=0; n < mesh._N; n++)
    {
        for(int i=0; i < rkParams.I; i++)
        {
            if( rkParams.lambda_i[i] != 0 )
                precursors[0][n][i] = rkParams.beta_i[i] * fissionProd[n]
                    / ( rkParams.lambda_i[i] * kcrit );
        }
    }

    // initialize new mesh
    Mesh newMesh = Mesh();

    // cycle through time steps
    for(int t=0; t < timeSteps.size()-1; t++)
    {
        // get current time step
        double time = timeSteps[t+1];
        double dt = timeSteps[t+1] - timeSteps[t];

        // interpolate mesh vector to form new mesh
        newMesh.interpolate(timeVector, meshVector, time);

        // check if new matrices need to be formed
        bool modified[] = {false, false, false, false};
        for(int n=0; n < mesh._N; n++)
        {
            // get materials
            XSdata * mat1 = mesh._material[n];
            XSdata * mat2 = newMesh._material[n];
            
            // check all groups
            for(int g=0; g < mesh._G; g++)
            {
                // check diffusion
                if( mat1->D[g] != mat2->D[g] )
                    modified[0] = true;
                
                // check absorption
                if( mat1->siga[g] != mat2->siga[g] )
                    modified[1] = true;

                // check scattering
                for(int gp=0; gp < mesh._G; gp++)
                    if( mat1->sigs[gp][g] != mat2->sigs[gp][g] )
                        modified[2] = true;

                // check fission
                if( mat1->nuSigf[g] != mat2->nuSigf[g] or 
                        mat1->chi[g] != mat2->chi[g])
                    modified[3] = true;
            }
        }

        // form new matrices as needed
        if(modified[0])
            Dhat = formDhatMatrix(newMesh, index);
        if(modified[1])
            SigA = formSigAMatrix(newMesh, index);
        if(modified[2])
            SigS = formSigSMatrix(newMesh, index);
        if(modified[3])
            F = formFMatrix(newMesh, index);

        // create tansient fission matrix
        Sparse Fhat = formFhatMatrix(newMesh, rkParams, dt, kcrit, index);

        // create T = (A - F) matrix
        Sparse T = Dhat + SigA + SigS - Fhat;

        // add time absorption term on diagonal
        for(int g=0; g < mesh._G; g++)
        {
            double time_abs = 1 / (rkParams.v[g] * dt);
            for(int n=0; n < mesh._N; n++)
            {
                double val = mesh._delta[n] * time_abs + 
                    T( index(n,g), index(n,g) );
                T.setVal( index(n,g), index(n,g), val );
            }
        }
        
        // copy newMesh to mesh
        mesh = newMesh;
        
        // create S vector
        std::vector<double> S = formSVector(mesh, rkParams, flux[t], 
                precursors[t], dt, kcrit, index);
        
        // solve flux
        flux[t+1] = T.optimalSOR(S, flux[t], tol, maxiters, sum_inner_iters);

        // tally total power (assuming nu is constant)
        double total_power = 0;

        // calculate new neutron precursors
        for(int n=0; n < mesh._N; n++)
        {
            // extract material
            XSdata* mat = mesh._material[n];

            // calculate total fission production
            double fission = 0;
            for(int g=0; g < mesh._G; g++)
                fission += mat->nuSigf[g] * flux[t+1][index(n,g)];

            total_power += fission;

            // calculate precursors for each group in I
            for(int i=0; i < rkParams.I; i++)
            {
                precursors[t+1][n][i] = (precursors[t][n][i] + 
                    rkParams.beta_i[i] * dt * fission / kcrit) /
                    (1 + rkParams.lambda_i[i] * dt);
            }
        }

        // add total power to power vector
        npower.push_back(total_power);
        std::cout << "Power = " << total_power << endl;
    }
    return;
}


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

    N = mesh._nodes;
   
    // form matrices
    /*
    Sparse A = Sparse(N*G, N*G);
    Sparse F = Sparse(N*G, N*G);
    formSteadyStateMatrixProblem(mesh, A, F, index);
    */
    Sparse F = formFMatrix(mesh, index);
    Sparse SigA = formSigAMatrix(mesh, index);
    Sparse SigS = formSigSMatrix(mesh, index);
    Sparse Dhat = formDhatMatrix(mesh, index);

    Sparse A = SigA + SigS + Dhat;

    /*
    Sparse C = A-F;
    std::cout << "A = " << endl;
    A.display();
    std::cout << "F = " << endl;
    F.display();
    std::cout << "C = " << endl;
    C.display();
    */

    // solve eigenvalue problem
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
            if(n == 0 or n == N-1)
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

// form Fhat transient matrix
Sparse formFhatMatrix(Mesh mesh, RKdata rkParams, double dt, double kcrit, 
        Indexer index)
{
    // extract dimensions
    int N = index._N;
    int G = index._G;

    // intialize matrix
    Sparse Fhat = Sparse(N*G, N*G);
    
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
                double prompt = (1 - rkParams.beta) / kcrit;

                double delayed = 0;
                for(int i=0; i<rkParams.I; i++)
                {
                    delayed += rkParams.beta_i[i] * rkParams.lambda_i[i] * dt
                        / ( (1 + rkParams.lambda_i[i] * dt) * kcrit );
                }
                
                double fission = (mat->chi[g] * prompt + 
                        rkParams.chi_d[g] * delayed) * mat->nuSigf[gp]*delta;
                
                if(fission != 0)
                    Fhat.setVal(index(n,g), index(n,gp), fission);
            }
        }
    }

    return Fhat;
}

// form S transient vector
std::vector<double> formSVector(Mesh mesh, RKdata rkParams, 
        std::vector<double> phi, std::vector<std::vector<double> > C, 
        double dt, double kcrit, Indexer index)
{
    // extract dimensions
    int N = index._N;
    int G = index._G;

    // intialize vector
    std::vector<double> S = std::vector<double>(N*G,0);
    
    // set vector terms
    for(int n=0; n<N; n++)
    {
        // calculate precursor contribution in cell
        double precursor = 0;
        for(int i=0; i < rkParams.I; i++)
        {
            precursor += rkParams.lambda_i[i] * C[n][i] 
                / (1 + rkParams.lambda_i[i] * dt);
        }

        for(int g=0; g<G; g++)
        {
            S[index(n,g)] = phi[index(n,g)] / (rkParams.v[g] * dt)
                    + rkParams.chi_d[g] * precursor;
        }
    }
    return S;
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
        
        // compute keff
        //k = sum_fission / (N*G);

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

/*
   Eigenvalue solutoin constructor
   */
eigenSolution::eigenSolution()
{
    keff = 0;
    outer_iters = 0;
    inner_iters = 0;
}

// destructor
eigenSolution::~eigenSolution() {};

/*
   creates 2D arrays of flux and power by group
   */
void eigenSolution::indexArrays(Indexer index)
{
    if( gFlux.size() > 0 )
    {
        std::cout << "Warning: arrays already indexed" << endl;
        return;
    }

    // fill arrays
    for(int g=0; g<index._G; g++)
    {
        std::vector<double> tempFlux;
        std::vector<double> tempPower;
        for(int n=0; n<index._N; n++)
        {
            tempFlux.push_back( flux[index(n,g)] );
            tempPower.push_back( power[index(n,g)] );
        }
        gFlux.push_back(tempFlux);
        gPower.push_back(tempPower);
    }
    return;
}

/*
   Eigenvalue solution destructor
   */       
std::vector<double> eigenSolution::getFlux(int g)
{
    if(gFlux.size() == 0)
    {
        std::cout << "Error: arrays not indexed!" << endl;
        return flux;
    }
    else
        return gFlux[g-1];
}

/*
   Function to return a vector of the fission reaction rate for group g
   */
std::vector<double> eigenSolution::getPower(int g)
{
    if(gPower.size() == 0)
    {
        std::cout << "Error: arrays not indexed!" << endl;
        return flux;
    }
    else
        return gPower[g-1];
} 

