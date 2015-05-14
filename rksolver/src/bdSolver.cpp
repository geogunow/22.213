#include"bdSolver.h"

//TODO: Description
rkSolution solveTransientBD(Transient trans, RKdata rkParams, int BD)
{
    //TODO: error checking
    std::cout << "Conducting BD Solver with BD = " << BD << endl;
     
    // get vectors
    std::vector<Mesh> meshVector = trans.meshVector;
    std::vector<double> timeVector = trans.timeVector;
    std::vector<double> timeSteps = trans.timeSteps;

    // set transient variables
    double tol = trans.tolerance;
    int maxiters = pow(10,6);
    double outer_tol = tol;
    double inner_tol = tol;
    int inner_solver = 2;
    int sum_inner_iters = 0;
   
    // get first mesh to solve steady state
    Mesh mesh = meshVector[0];

    // initialie matrices
    Indexer index = Indexer(mesh._N, mesh._G);

    // set N to the number of nodes
    int N = mesh._nodes;
 
    // intilaize solution vecotrs
    std::vector<std::vector<std::vector<double> > > precursors;
    std::vector<std::vector<double> > flux;
    for(int b=0; b < BD; b++)
    {
        std::vector<double> temp_flux(mesh._N * mesh._G, 0);
        std::vector<std::vector<double> > temp_precursors;
        for(int n=0; n < mesh._N; n++)
        {
            std::vector<double> temp_precursors2(rkParams.I, 0);
            temp_precursors.push_back(temp_precursors2);
        }
        flux.push_back(temp_flux);
        precursors.push_back(temp_precursors);
    }

    // initialize vectors to store weighted flux and precursor sums
    std::vector<double> F_sum(mesh._N * mesh._G, 0);
    std::vector<std::vector<double> > C_sum;
    for(int n=0; n < mesh._N; n++)
    {
        std::vector<double> temp(rkParams.I, 0);
        C_sum.push_back(temp);
    }
    std::vector<double> F_temp;
    std::vector<std::vector<double> > C_temp;

    // create transient reactor kinetics soution object
    rkSolution rkResult = rkSolution();

    // Assign BD alphas and omegas
    std::vector<double> alpha(BD+1);
    double omega;
    switch(BD){
        case 1:
        {
            double div = 1;
            double alpha_vals[] = {div, -1};
            for(int b=0; b < BD+1; b++)
                alpha[b] = alpha_vals[b] / div;
            omega = 1;
            break;
        }
        case 2:
        {
            double div = 3;
            double alpha_vals[] = {div, -4, 1};
            for(int b=0; b < BD+1; b++)
                alpha[b] = alpha_vals[b] / div;
            omega = 2.0 / div;
            break;
        }
        case 3:
        {
            double div = 11;
            double alpha_vals[] = {div, -18, 9, -2};
            for(int b=0; b < BD+1; b++)
                alpha[b] = alpha_vals[b] / div;
            omega = 6.0 / div;
            break;
        }
        case 4:
        {
            double div = 25;
            double alpha_vals[] = {div, -48, 36, -16, 3};
            for(int b=0; b < BD+1; b++)
                alpha[b] = alpha_vals[b] / div;
            omega = 12.0 / div;
            break;
        }
        case 5:
        {
            double div = 137;
            double alpha_vals[] = {div, -300, 300, -200, 75, -12};
            for(int b=0; b < BD+1; b++)
                alpha[b] = alpha_vals[b] / div;
            omega = 60.0 / div;
            break;
        }
        case 6:
        {
            double div = 147;
            double alpha_vals[] = {div, -360, 450, -400, 225, -72, 10};
            for(int b=0; b < BD+1; b++)
                alpha[b] = alpha_vals[b] / div;
            omega = 60.0 / div;
            break;
        }
        default:
            std::cout << "Error: BD order not supported" << endl;
            std::cout << "Please choose an order between 1 and 6" << endl;
            return rkResult;
    }

   
    // form steady state matrices
    Sparse F = formFMatrix(mesh, index);
    Sparse SigA = formSigAMatrix(mesh, index);
    Sparse SigS = formSigSMatrix(mesh, index);
    Sparse Dhat = formDhatMatrix(mesh, index);
    Sparse A = SigA + SigS + Dhat;
    
    // solve eigenvalue problem
    eigenSolution ssResult = eigen_solver(A, F, N, mesh._G, outer_tol, 
            inner_tol, maxiters, inner_solver);
    
    // extract critical k for transient solution
    double kcrit = ssResult.keff;

    // initialize first time step of flux
    flux[0] = ssResult.flux;

    // total fission tally
    double tot_fission = 0;
    
    // calculate intial precursor population for each cell
    for(int n=0; n < mesh._N; n++)
    {
        // calculate fission production in the cell
        double fp = 0;
        for(int g=0; g < mesh._G; g++)
            fp += ssResult.flux[index(n,g)] * mesh._material[n]->nuSigf[g]; 

        // add to total fission tally
        tot_fission += fp;

        // calculate initial precursor population
        for(int i=0; i < rkParams.I; i++)
        {
            if( rkParams.lambda_i[i] != 0 )
                precursors[0][n][i] = rkParams.beta_i[i] * fp
                    / ( rkParams.lambda_i[i] * kcrit );
        }
    }

    // record power and profile
    rkResult.power.push_back(tot_fission);
    rkResult.powerProfile.push_back(ssResult.power);

    // cycle through time steps
    for(int t=0; t < timeSteps.size()-1; t++)
    {
        // get current time step
        double time = timeSteps[t+1];
        double dt = timeSteps[t+1] - timeSteps[t];

        // initialize new mesh
        Mesh newMesh = Mesh();
        
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
                if( mat1->nuSigf[g] != mat2->nuSigf[g] || 
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
        Sparse Fhat = formFhatMatrixBD(newMesh, rkParams, dt, kcrit, index,
                alpha[0], omega);

        // create T = (A - F) matrix
        Sparse T = Dhat + SigA + SigS - Fhat;

        // add time absorption term on diagonal
        for(int g=0; g < mesh._G; g++)
        {
            double time_abs = alpha[0] / (omega * rkParams.v[g] * dt);
            for(int n=0; n < mesh._N; n++)
            {
                double val = mesh._delta[n] * time_abs + 
                    T( index(n,g), index(n,g) );
                T.setVal( index(n,g), index(n,g), val );
            }
        }

        // calculate weighted sum of precursors and flux
        for(int n=0; n < mesh._N; n++)
        {
            for(int g=0; g < mesh._G; g++)
            {
                // calculate sum of flux accross back vectors
                F_sum[index(n,g)] = 0;
                for(int b=0; b < BD; b++)
                {
                    // check for beginning cases not to use back vectors
                    if(t < b)
                        F_temp = flux[0];
                    else
                        F_temp = flux[b];

                    F_sum[index(n,g)] += alpha[b+1] * F_temp[index(n,g)];
                }
            }

            // repeat for precursors
            for(int i=0; i < rkParams.I; i++)
            {
                // calculate sum of precursors accross back vectors
                C_sum[n][i] = 0;
                for(int b=0; b < BD; b++)
                {
                    // check for back vectors not available
                    if(t < b)
                        C_temp = precursors[0];
                    else
                        C_temp = precursors[b];
                    C_sum[n][i] += alpha[b+1] * C_temp[n][i];
                }
            }
        }
                
        // create S vector
        std::vector<double> S = formSVectorBD(mesh, rkParams, F_sum, C_sum, dt,
                kcrit, index, alpha[0], omega);
        
        // copy newMesh to mesh
        mesh = newMesh;
        
        // copy flux and precursors to back vectors
        for(int b=BD-1; b > 0; b--)
        {
            flux[b] = flux[b-1];
            precursors[b] = precursors[b-1];
        }

        // solve flux
        flux[0] = T.optimalSOR(S, flux[0], tol, maxiters, sum_inner_iters);

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
                fission += mat->nuSigf[g] * flux[0][index(n,g)];

            total_power += fission;

            // calculate precursors for each group in I
            for(int i=0; i < rkParams.I; i++)
            {
                double lambda_i = rkParams.lambda_i[i];
                double beta_i = rkParams.beta_i[i];
                precursors[0][n][i] = (beta_i * omega * dt * fission / kcrit 
                        - C_sum[n][i]) / (alpha[0] + lambda_i * omega * dt);
            }
        }

        // add total power to power vector
        rkResult.power.push_back(total_power);
        
        if( (t+2)%100 == 0 or (t+2) == timeSteps.size())
            std::cout << "Completed " << t+2 << "/" << timeSteps.size()
                << " timesteps" << endl;
    }
    return rkResult;
}

// form Fhat transient matrix
Sparse formFhatMatrixBD(Mesh mesh, RKdata rkParams, double dt, double kcrit, 
        Indexer index, double alpha, double omega)
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
                    double lambda_i = rkParams.lambda_i[i];
                    double beta_i = rkParams.beta_i[i];
                    delayed += beta_i * lambda_i * omega * dt 
                        / (kcrit * (alpha + lambda_i * omega * dt));
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
std::vector<double> formSVectorBD(Mesh mesh, RKdata rkParams, 
        std::vector<double> F_sum, std::vector<std::vector<double> > C_sum,
        double dt, double kcrit, Indexer index, double alpha, double omega)
{
    // extract dimensions
    int N = index._N;
    int G = index._G;

    // intialize vector
    std::vector<double> S = std::vector<double>(N*G,0);
    
    // set vector terms
    for(int n=0; n < N; n++)
    {
        // get mesh spacing
        double dx = mesh._delta[n];

        // calculate precursor contribution in cell
        double precursor = 0;
        for(int i=0; i < rkParams.I; i++)
        {
            double lambda_i = rkParams.lambda_i[i];
            precursor += lambda_i * C_sum[n][i] 
                / (alpha + lambda_i * omega * dt);
        }

        for(int g=0; g<G; g++)
        {
            S[index(n,g)] = -1 * (F_sum[index(n,g)] * dx 
                    / (rkParams.v[g] * omega * dt)
                    + rkParams.chi_d[g] * precursor * dx);
        }
    }
    return S;
}

