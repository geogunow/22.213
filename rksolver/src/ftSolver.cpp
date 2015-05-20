#include"ftSolver.h"

//TODO: Description
rkSolution solveTransientFT(Transient trans, RKdata rkParams, int omegaMode)
{
    // get vectors
    std::vector<Mesh> meshVector = trans.meshVector;
    std::vector<double> timeVector = trans.timeVector;
    std::vector<double> timeSteps = trans.timeSteps;

    // set transient variables
    double tol = trans.tolerance;
    int maxiters = pow(10,6);
    double outer_tol = tol;
    double inner_tol = tol;
    int inner_solver = 1;
    int sum_inner_iters = 0;

    //TODO: error checking
    
    // create transient reactor kinetics soution object
    rkSolution rkResult = rkSolution();

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
    eigenSolution ssResult = eigen_solver(A, F, N, mesh._G, outer_tol, 
            inner_tol, maxiters, inner_solver);
   
    // extract critical k for transient solution
    double kcrit = ssResult.keff;

    // initialize omega to be zeros
    std::vector<double> omega(mesh._N, 0);

    // intialize precursor and flux structure
    std::vector<std::vector<double> > precursors;
    std::vector<double> flux(mesh._N * mesh._G, 0);
    for(int n=0; n < mesh._N; n++)
    {
        std::vector<double> temp (rkParams.I, 0);
        precursors.push_back(temp);
    }
    
    // initialize first time step of flux
    flux = ssResult.flux;

    // Normalize intial flux
    double tot_flux = 0;
    for(int i=0; i < flux.size(); i++)
        tot_flux += ssResult.flux[i];
    for(int i=0; i < flux.size(); i++)
        flux[i] /= tot_flux;

    // total fission tally
    double tot_fission = 0;
    
    // calculate intial precursor population for each cell
    for(int n=0; n < mesh._N; n++)
    {
        // calculate fission production in the cell
        double fp = 0;
        for(int g=0; g < mesh._G; g++)
            fp += flux[index(n,g)] * mesh._material[n]->nuSigf[g]; 

        // add to total fission tally
        tot_fission += fp * mesh._delta[n];

        // calculate initial precursor population
        for(int i=0; i < rkParams.I; i++)
        {
            if( rkParams.lambda_i[i] != 0 )
                precursors[n][i] = rkParams.beta_i[i] * fp
                    / ( rkParams.lambda_i[i] * kcrit );
        }
    }

    // boolean variable for synchronized mode
    bool firstIter = true;

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
        Sparse Fhat = formFhatMatrixFT(F, newMesh, rkParams, dt, kcrit, 
                omega, index);

        // create T = (A - F) matrix
        Sparse T = Dhat + SigA + SigS - Fhat;

        // add time absorption term on diagonal
        for(int g=0; g < mesh._G; g++)
        {
            double time_abs = 1 / (rkParams.v[g] * dt);
            for(int n=0; n < mesh._N; n++)
            {
                double val = mesh._delta[n] * time_abs * (1 + omega[n]*dt) + 
                    T( index(n,g), index(n,g) );
                T.setVal( index(n,g), index(n,g), val );
            }
        }
       
        // create S vector
        std::vector<double> S = formSVectorFT(newMesh, rkParams, F*flux, flux, 
                precursors, dt, kcrit, omega, index);
 
        // copy newMesh to mesh
        mesh = newMesh;
        
        // solve flux
        sum_inner_iters = 0;
        std::vector<double> psi = 
            T.gaussSeidel(S, flux, tol, maxiters, sum_inner_iters);

        // record inner iterations
        rkResult.innerSteps.push_back(sum_inner_iters);

        // calculate old and new total fluxes and reassign fluxes
        std::vector<double> old_flux_total(mesh._N,0);
        std::vector<double> new_flux_total(mesh._N,0);
        for(int n=0; n < mesh._N; n++)
        {
            for(int g=0; g < mesh._G; g++)
            {
                old_flux_total[n] += flux[index(n,g)];
                new_flux_total[n] += psi[index(n,g)] * exp(omega[n]*dt);
            }
        }

        // update omega globally if requested
        if(omegaMode == 1 or omegaMode == 3)
        {
            double sum_old_flux = 0, sum_new_flux = 0;
            for(int n=0; n < mesh._N; n++)
            {
                sum_old_flux += old_flux_total[n];
                sum_new_flux += new_flux_total[n];
            }
            for(int n=0; n < mesh._N; n++)
                omega[n] = log(sum_new_flux / sum_old_flux) / dt;
        }

        // update omega locally if requested
        else if(omegaMode == 2 or omegaMode == 4)
            for(int n=0; n < mesh._N; n++)
                omega[n] = log(new_flux_total[n] / old_flux_total[n]) / dt;

        // copy new flux estimate if unsynchronized
        if(omegaMode == 1 or omegaMode == 2
                or (omegaMode > 2 && !firstIter) )
        {
            // calculate new neutron precursors
            for(int n=0; n < mesh._N; n++)
            {
                // extract material
                XSdata* mat = mesh._material[n];

                // calculate for each delayed group
                for(int i=0; i < rkParams.I; i++)
                {
                    // load decay constant
                    double lambda_i = rkParams.lambda_i[i];

                    // calculate decayed precursors
                    precursors[n][i] = precursors[n][i] * exp(-lambda_i * dt);

                    // add production of precursors
                    for(int g=0; g < mesh._G; g++)
                    {
                        double Ag = rkParams.beta_i[i] * mat->nuSigf[g] /
                            (kcrit * dt * (lambda_i + omega[n]));
                        double Bg = (psi[index(n,g)] - flux[index(n,g)]) /
                            (lambda_i + omega[n]);

                        precursors[n][i] += Ag * ( (psi[index(n,g)] * dt - Bg)
                                * exp(omega[n] * dt) - (flux[index(n,g)] * dt - Bg)
                                * exp(-lambda_i * dt) );
                    }
                }
            }
            // copy flux
            for(int n=0; n < mesh._N; n++)
                for(int g=0; g < mesh._G; g++)
                   flux[index(n,g)] = psi[index(n,g)] * exp(omega[n]*dt);
        }

        // calculate total power (assuming nu is constant)
        double total_power = 0;
        std::vector<double> power = F*flux;
        for(int n=0; n < mesh._N; n++)
            for(int g=0; g < mesh._G; g++)
                total_power += power[index(n,g)];
       
        // roll back solution if synchronized mode requested
        if(omegaMode > 2 and firstIter)
        {
            // roll back time step
            t -= 1;
            firstIter = false;
        }
        else
        {
            // add total power to power vector
            rkResult.power.push_back(total_power);
            firstIter = true;
        }
        
        if( (t+2)%100 == 0 or (t+2) == timeSteps.size())
            std::cout << "Completed " << t+2 << "/" << timeSteps.size()
                << " timesteps" << endl;
    }
    return rkResult;
}

// form Fhat frequency transform transient matrix
Sparse formFhatMatrixFT(Sparse F, Mesh mesh, RKdata rkParams, double dt, 
        double kcrit, std::vector<double> omega, Indexer index)
{
    // extract dimensions
    int N = index._N;
    int G = index._G;

    // intialize matrix
    Sparse Fhat = Sparse(N*G, N*G);
    
    // setup matrix
    for(int n=0; n<N; n++)
    {
        XSdata* mat = mesh._material[n];
        for(int g=0; g<G; g++)
        {
            // extract spectrum
            double chi = mat->chi[g];

            // calculate fission matrix terms
            for(int gp=0; gp<G; gp++)
            {
                if(F(index(n,g), index(n,gp)) != 0)
                {
                    double prompt = (1 - rkParams.beta) * mat->chi[g] 
                        / mat->chi[g];

                    double delayed = 0;
                    for(int i=0; i<rkParams.I; i++)
                    {
                        double lambda_i = rkParams.lambda_i[i];
                        double beta_i = rkParams.beta_i[i];
                        double chi_d = rkParams.chi_d[g];
                        delayed += chi_d * beta_i * lambda_i / (chi
                            * dt * (lambda_i + omega[n])) * (dt 
                            - 1.0 / (lambda_i + omega[n]) 
                            + exp(-(lambda_i + omega[n]) * dt) 
                            / (lambda_i + omega[n]));
                    }
                    
                    double fission = F(index(n,g), index(n,gp)) 
                        * (prompt + delayed) / kcrit;
                    Fhat.setVal(index(n,g), index(n,gp), fission);
                }
            }
        }
    }

    return Fhat;
}

//TODO: wrtie description
std::vector<double> formSVectorFT(Mesh mesh, RKdata rkParams, 
        std::vector<double> source, std::vector<double> flux, 
        std::vector<std::vector<double> > C, double dt, double kcrit, 
        std::vector<double> omega, Indexer index)
{
    // extract dimensions
    int N = index._N;
    int G = index._G;

    // intialize vector
    std::vector<double> S = std::vector<double>(N*G,0);
    
    // set vector terms
    for(int n=0; n<N; n++)
    {
        // get mesh spacing
        double dx = mesh._delta[n];

        // loop over groups
        for(int g=0; g<G; g++)
        {
            // get spectrum
            double chi = mesh._material[n]->chi[g];
        
            // calculate precursor contribution in cell
            double precursor = 0;
            
            for(int i=0; i < rkParams.I; i++)
            {
                double chi_d = rkParams.chi_d[g];
                double lambda_i = rkParams.lambda_i[i];
                double beta_i = rkParams.beta_i[i];
                
                precursor += chi_d * lambda_i * C[n][i] * 
                        exp(-(lambda_i + omega[n]) * dt);
                if( chi != 0 )
                {
                    double q1 = chi_d * lambda_i * beta_i * source[index(n,g)]
                        / (chi * dx * kcrit * (lambda_i + omega[n]));
                    double experm = exp(-(lambda_i + omega[n]) * dt);
                    precursor += q1 * ( -experm + (1-experm)
                            / ((lambda_i + omega[n]) * dt) );
                }
            }

            S[index(n,g)] = flux[index(n,g)] * dx / (rkParams.v[g] * dt)
                    + precursor * dx;
        }
    }
    return S;
}
