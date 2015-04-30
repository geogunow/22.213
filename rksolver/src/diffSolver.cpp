#include"diffSolver.h"

//TODO: Description
rkSolution solveTransient(Transient trans, RKdata rkParams)
{
    // get vectors
    std::vector<Mesh> meshVector = trans.meshVector;
    std::vector<double> timeVector = trans.timeVector;
    std::vector<double> timeSteps = trans.timeSteps;

    // set transient variables
    double tol = pow(10,-6);
    int maxiters = pow(10,6);
    double outer_tol = pow(10,-6);
    double inner_tol = pow(10,-6);
    int inner_solver = 2;
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
                precursors[n][i] = rkParams.beta_i[i] * fp
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
        
        // create S vector
        std::vector<double> S = formSVector(mesh, rkParams, flux, 
                precursors, dt, kcrit, index);
        
        // copy newMesh to mesh
        mesh = newMesh;


        // solve flux
        flux = T.optimalSOR(S, flux, tol, maxiters, sum_inner_iters);

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
                fission += mat->nuSigf[g] * flux[index(n,g)];

            total_power += fission;

            // calculate precursors for each group in I
            for(int i=0; i < rkParams.I; i++)
            {
                precursors[n][i] = (precursors[n][i] + 
                    rkParams.beta_i[i] * dt * fission / kcrit) /
                    (1 + rkParams.lambda_i[i] * dt);
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


//TODO: Description
//TODO: Make better adjoint weighting input
rkSolution solvePKE(Transient trans, RKdata rkParams, int shape_step,
        int adj_weighting, int adj_step)
{
    // get vectors
    std::vector<Mesh> meshVector = trans.meshVector;
    std::vector<double> timeVector = trans.timeVector;
    std::vector<double> timeSteps = trans.timeSteps;

    // set transient variables
    double tol = pow(10,-6);
    int maxiters = pow(10,6);
    double outer_tol = pow(10,-8);
    double inner_tol = pow(10,-8);
    int inner_solver = 2;
    int sum_inner_iters = 0;

    //TODO: error checking
    if(shape_step >= meshVector.size())
        std::cout << "Error: Shape step exceeds mesh vector size" << endl;
    
    if(adj_step >= meshVector.size())
        std::cout << "Error: Adjoint step exceeds mesh vector size" << endl;
    
    // create transient reactor kinetics soution object
    rkSolution rkResult = rkSolution();

    // get the mesh for the shape time step
    Mesh mesh = meshVector[shape_step];
    
    // initialize matrices
    Indexer index = Indexer(mesh._N, mesh._G);

    // set N to the number of nodes
    int N = mesh._nodes;
   
    // solve eigenvalue problem
    eigenSolution ssResult = solveCritical(mesh, outer_tol, inner_tol, 
            maxiters, inner_solver);

    // extract shape function for transient solution
    std::vector<double> shape = ssResult.flux;

    // extract critical k for transient solution
    eigenSolution result0 = solveCritical(meshVector[0], outer_tol, inner_tol,
            maxiters, inner_solver);
    double kcrit = result0.keff;

    // extract shape power from shape function
    std::vector<double> shape_power (mesh._G, 0);
    for(int g=0; g < mesh._G; g++)
        for(int n=0; n < mesh._N; n++)
            shape_power[g] += ssResult.power[index(n,g)];
   
    // get the function for adjoint weighting
    std::vector<double> adjoint (shape.size(), 1);
    
    // solve adjoint problem if necessary
    if(adj_weighting == 1)
    {
        int L = mesh._N * mesh._G;
        Sparse Aadj = Sparse(L,L);
        Sparse Fadj = Sparse(L,L);
        formSteadyStateMatrixProblem(meshVector[adj_step], Aadj, Fadj, index);
        Sparse At = Aadj.transpose();
        Sparse Ft = Fadj.transpose();
        eigenSolution adjResult = eigen_solver(At, Ft, N, mesh._G, 
                outer_tol, inner_tol, maxiters, inner_solver);
        adjoint = adjResult.flux;
    }
    
    // intialize precursor and flux structure
    std::vector<std::vector<double> > C_tilde;
    std::vector<double> flux(mesh._G, 0);
    for(int g=0; g < mesh._G; g++)
    {
        std::vector<double> temp (rkParams.I, 0);
        C_tilde.push_back(temp);
    }
    
    // total fission tally
    double tot_fission = 0;
    
    // calculate initial precursors
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
            {
                double precursors = rkParams.beta_i[i] * fp
                    / ( rkParams.lambda_i[i] * kcrit );
                
                for(int g=0; g < mesh._G; g++)
                    C_tilde[g][i] += rkParams.chi_d[g] * mesh._delta[n]
                        * adjoint[index(n,g)] * precursors;
            }
        }
    }
    
    // precalculate time absorption without dt
    std::vector<double> time_abs(mesh._G, 0);
    for(int g=0; g < mesh._G; g++)
        for(int n=0; n < mesh._N; n++)
            time_abs[g] += mesh._delta[n] * adjoint[index(n,g)]
                * shape[index(n,g)] / rkParams.v[g];
    
    // demand unity initial amplitudes
    std::vector<double> amplitude (mesh._G, 1);

    // form PKE matricies
    Sparse SigS = formSigSMatrixPKE(mesh, index, shape, adjoint);
    Sparse SigA = formSigAMatrixPKE(mesh, index, shape, adjoint);
    Sparse Dhat = formDhatMatrixPKE(mesh, index, shape, adjoint);
    Sparse F = formFMatrixPKE(mesh, index, shape, adjoint, rkParams, kcrit, 
            timeSteps[1] - timeSteps[0]);
    
    // setup matrix for time dependent problem
    Sparse T = SigA + SigS + Dhat - F;
    
    // add time absorption term on diagonal
    double dt = timeSteps[1] - timeSteps[0];
    for(int g=0; g < mesh._G; g++)
    {
        double val = time_abs[g] / dt + T(g,g);
        T.setVal(g,g,val);
    }

    // record power and profile
    rkResult.powerProfile.push_back(amplitude);
    rkResult.power.push_back( dot(shape_power, amplitude) );

    // cycle through time steps
    for(int t=0; t < timeSteps.size()-1; t++)
    {
        // get current time step
        double time = timeSteps[t+1];
        dt = timeSteps[t+1] - timeSteps[t];

        // initialize new mesh
        Mesh newMesh = Mesh();
        
        // interpolate mesh vector to form new mesh
        newMesh.interpolate(timeVector, meshVector, time);

        // check if new matrices need to be formed
        bool modified[] = {false, false, false, false,false};
        for(int n=0; n < mesh._N; n++)
        {
            // get materials
            XSdata * mat1 = mesh._material[n];
            XSdata * mat2 = newMesh._material[n];
            
            // check all groups
            for(int g=0; g < mesh._G; g++)
            {
                // check absorption
                if( mat1->siga[g] != mat2->siga[g] )
                {
                    modified[0] = true;
                    modified[1] = true;
                }
                // check scattering
                for(int gp=0; gp < mesh._G; gp++)
                    if( mat1->sigs[gp][g] != mat2->sigs[gp][g] )
                    {
                        modified[0] = true;
                        modified[2] = true;
                    }
                // check fission
                if( mat1->nuSigf[g] != mat2->nuSigf[g] || 
                        mat1->chi[g] != mat2->chi[g])
                {
                    modified[0] = true;
                    modified[3] = true;
                }
                // check diffusion
                if( mat1->D[g] != mat2->D[g])
                {
                    modified[0] = true;
                    modified[4] = true;
                }
            }
        }

        // test for changed time step
        double time_step_change;
        if(t > 0)
        {
            double time_step1 = timeSteps[t+1] - timeSteps[t];
            double time_step2 = timeSteps[t] - timeSteps[t-1];
            time_step_change = abs(time_step1 - time_step2)/time_step2;
        }
        else
            time_step_change = 0;

        // form new matrices as needed
        if(modified[0] or time_step_change > pow(10,-6))
        {
            if(modified[1])
                SigA = formSigAMatrixPKE(newMesh, index, shape, adjoint);
        
            if(modified[2])
                SigS = formSigSMatrixPKE(newMesh, index, shape, adjoint);
        
            if(modified[3] or time_step_change > pow(10,-6))
                F = formFMatrixPKE(newMesh, index, shape, adjoint, 
                        rkParams, kcrit, dt);
            
            if(modified[4])
                Dhat = formDhatMatrixPKE(newMesh, index, shape, adjoint);
        
            // recreate T = (A - F) matrix
            T = SigA + SigS + Dhat - F;

            // add time absorption term on diagonal
            for(int g=0; g < mesh._G; g++)
            {
                double val = time_abs[g] / dt + T(g,g);
                T.setVal(g,g,val);
            }
        }

        // create S vector
        std::vector<double> S = formSVectorPKE(index, amplitude, rkParams, 
                C_tilde, dt, time_abs);
    
        // copy newMesh to mesh
        mesh = newMesh;
        
        // solve flux
        amplitude = T.optimalSOR(S, amplitude, tol, maxiters, sum_inner_iters);

        // precalcualte R_pi term
        std::vector<double> rpi(mesh._G, 0);
        for(int g=0; g < mesh._G; g++)
            for(int gp=0; gp < mesh._G; gp++)
                for(int n=0; n < mesh._N; n++)
                {
                    XSdata* mat = mesh._material[n];
                    rpi[g] += rkParams.chi_d[g] * adjoint[index(n,g)]
                        * mat->nuSigf[gp] * shape[index(n,gp)] * amplitude[gp]
                        * mesh._delta[n];
                }

        // calculate new precursors
        for(int g=0; g < mesh._G; g++)
        {
            for(int i=0; i < rkParams.I; i++)
            {
                C_tilde[g][i] += rkParams.beta_i[i] * dt / kcrit * rpi[g];
                C_tilde[g][i] /= (1 + rkParams.lambda_i[i] * dt);
            }
        }

        // add total power to power vector
        rkResult.powerProfile.push_back(amplitude);
        rkResult.power.push_back( dot(shape_power, amplitude) );

        if( (t+2)%100 == 0 or (t+2) == timeSteps.size())
            std::cout << "Completed " << t+2 << "/" << timeSteps.size()
                << " timesteps" << endl;
    }
    return rkResult;
}

//TODO: Description
//TODO: Make better adjoint weighting input
rkSolution solvePKESimple(Transient trans, RKdata rkParams, int shape_step,
        int adj_weighting, int adj_step)
{
    // get vectors
    std::vector<Mesh> meshVector = trans.meshVector;
    std::vector<double> timeVector = trans.timeVector;
    std::vector<double> timeSteps = trans.timeSteps;

    // set transient variables
    double tol = pow(10,-6);
    int maxiters = pow(10,6);
    double outer_tol = pow(10,-8);
    double inner_tol = pow(10,-8);
    int inner_solver = 2;
    int sum_inner_iters = 0;

    //TODO: error checking
    if(shape_step >= meshVector.size())
        std::cout << "Error: Shape step exceeds mesh vector size" << endl;
    
    if(adj_step >= meshVector.size())
        std::cout << "Error: Adjoint step exceeds mesh vector size" << endl;
    
    // create transient reactor kinetics soution object
    rkSolution rkResult = rkSolution();

    // get the mesh for the shape time step
    Mesh mesh = meshVector[shape_step];
    
    // initialize matrices
    Indexer index = Indexer(mesh._N, mesh._G);

    // set N to the number of nodes
    int N = mesh._nodes;
    
    // solve eigenvalue problem
    eigenSolution ssResult = solveCritical(mesh, outer_tol, inner_tol, 
            maxiters, inner_solver);

    // extract shape function for transient solution
    std::vector<double> shape = ssResult.flux;

    // extract critical k for transient solution
    eigenSolution result0 = solveCritical(meshVector[0], outer_tol, inner_tol,
            maxiters, inner_solver);
    double kcrit = result0.keff;

    // extract shape power from shape function
    double shape_power = 0;
    for(int g=0; g < mesh._G; g++)
        for(int n=0; n < mesh._N; n++)
            shape_power += ssResult.power[index(n,g)];
   
    // get the function for adjoint weighting
    std::vector<double> adjoint (shape.size(), 1);
    
    // solve adjoint problem if necessary
    if(adj_weighting == 1)
    {
        int L = mesh._N * mesh._G;
        Sparse Aadj = Sparse(L,L);
        Sparse Fadj = Sparse(L,L);
        formSteadyStateMatrixProblem(meshVector[adj_step], Aadj, Fadj, index);
        Sparse At = Aadj.transpose();
        Sparse Ft = Fadj.transpose();
        eigenSolution adjResult = eigen_solver(At, Ft, N, mesh._G, 
                outer_tol, inner_tol, maxiters, inner_solver);
        adjoint = adjResult.flux;
    }
    
    // intialize precursor and flux structure
    std::vector<double> C_tilde(rkParams.I,0);
    
    // total fission tally
    double tot_fission = 0;
    
    // calculate initial precursors
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
            {
                double precursors = rkParams.beta_i[i] * fp
                    / ( rkParams.lambda_i[i] * kcrit );
                
                for(int g=0; g < mesh._G; g++)
                    C_tilde[i] += rkParams.chi_d[g] * mesh._delta[n]
                        * adjoint[index(n,g)] * precursors;
            }
        }
    }
    
    
    // precalculate time absorption without dt
    double time_abs = 0;
    for(int g=0; g < mesh._G; g++)
        for(int n=0; n < mesh._N; n++)
            time_abs += mesh._delta[n] * adjoint[index(n,g)]
                * shape[index(n,g)] / rkParams.v[g];
    
    // demand unity initial amplitudes
    double amplitude = 1;

    // form PKE matricies
    Sparse SigS = formSigSMatrixPKE(mesh, index, shape, adjoint);
    Sparse SigA = formSigAMatrixPKE(mesh, index, shape, adjoint);
    Sparse Dhat = formDhatMatrixPKE(mesh, index, shape, adjoint);
    Sparse F = formFMatrixPKE(mesh, index, shape, adjoint, rkParams, kcrit, 
            timeSteps[1] - timeSteps[0]);
    
    // setup matrix for time dependent problem
    Sparse T = SigA + SigS + Dhat - F;
    
    // calcualte time step
    double dt = timeSteps[1] - timeSteps[0];
    
    // collapse matrix
    double lhs = T.sum() + time_abs / dt;

    // record power and profile
    rkResult.power.push_back(shape_power * amplitude);

    // cycle through time steps
    for(int t=0; t < timeSteps.size()-1; t++)
    {
        // get current time step
        double time = timeSteps[t+1];
        dt = timeSteps[t+1] - timeSteps[t];

        // initialize new mesh
        Mesh newMesh = Mesh();
        
        // interpolate mesh vector to form new mesh
        newMesh.interpolate(timeVector, meshVector, time);

        // check if new matrices need to be formed
        bool modified[] = {false, false, false, false,false};
        for(int n=0; n < mesh._N; n++)
        {
            // get materials
            XSdata * mat1 = mesh._material[n];
            XSdata * mat2 = newMesh._material[n];
            
            // check all groups
            for(int g=0; g < mesh._G; g++)
            {
                // check absorption
                if( mat1->siga[g] != mat2->siga[g] )
                {
                    modified[0] = true;
                    modified[1] = true;
                }
                // check scattering
                for(int gp=0; gp < mesh._G; gp++)
                    if( mat1->sigs[gp][g] != mat2->sigs[gp][g] )
                    {
                        modified[0] = true;
                        modified[2] = true;
                    }
                // check fission
                if( mat1->nuSigf[g] != mat2->nuSigf[g] || 
                        mat1->chi[g] != mat2->chi[g])
                {
                    modified[0] = true;
                    modified[3] = true;
                }
                // check diffusion
                if( mat1->D[g] != mat2->D[g])
                {
                    modified[0] = true;
                    modified[4] = true;
                }
            }
        }

        // test for changed time step
        double time_step_change;
        if(t > 0)
        {
            double time_step1 = timeSteps[t+1] - timeSteps[t];
            double time_step2 = timeSteps[t] - timeSteps[t-1];
            time_step_change = abs(time_step1 - time_step2)/time_step2;
        }
        else
            time_step_change = 0;

        // form new matrices as needed
        if(modified[0] or time_step_change > pow(10,-6))
        {
            if(modified[1])
                SigA = formSigAMatrixPKE(newMesh, index, shape, adjoint);
        
            if(modified[2])
                SigS = formSigSMatrixPKE(newMesh, index, shape, adjoint);
        
            if(modified[3] or time_step_change > pow(10,-6))
                F = formFMatrixPKE(newMesh, index, shape, adjoint, 
                        rkParams, kcrit, dt);
            
            if(modified[4])
                Dhat = formDhatMatrixPKE(newMesh, index, shape, adjoint);
        
            // recreate T = (A - F) matrix
            T = SigA + SigS + Dhat - F;

            // collapse matrix
            lhs = T.sum() + time_abs / dt;
        }

        // form right hand side
        double rhs = amplitude * time_abs / dt;
        for(int i=0; i < rkParams.I; i++)
            rhs += rkParams.lambda_i[i] * C_tilde[i]
                / (1 + rkParams.lambda_i[i] * dt);

        // copy newMesh to mesh
        mesh = newMesh;
        
        // solve flux
        amplitude = rhs / lhs;

        // precalcualte R_pi term
        double rpi = 0;
        for(int g=0; g < mesh._G; g++)
            for(int gp=0; gp < mesh._G; gp++)
                for(int n=0; n < mesh._N; n++)
                {
                    XSdata* mat = mesh._material[n];
                    rpi += rkParams.chi_d[g] * adjoint[index(n,g)]
                        * mat->nuSigf[gp] * shape[index(n,gp)] * amplitude
                        * mesh._delta[n];
                }

        // calculate new precursors
        for(int i=0; i < rkParams.I; i++)
        {
            C_tilde[i] += rkParams.beta_i[i] * dt / kcrit * rpi;
            C_tilde[i] /= (1 + rkParams.lambda_i[i] * dt);
        }

        // add total power to power vector
        rkResult.power.push_back(shape_power * amplitude);

        if( (t+2)%100 == 0 or (t+2) == timeSteps.size())
            std::cout << "Completed " << t+2 << "/" << timeSteps.size()
                << " timesteps" << endl;
    }
    return rkResult;
}


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
                psi[index(n,g)] *= exp(omega[n]*dt);
                new_flux_total[n] += psi[index(n,g)];
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
            // copy flux
            for(int n=0; n < mesh._N; n++)
                for(int g=0; g < mesh._G; g++)
                   flux[index(n,g)] = psi[index(n,g)];
 
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
        // get mesh spacing
        double dx = mesh._delta[n];

        // calculate precursor contribution in cell
        double precursor = 0;
        for(int i=0; i < rkParams.I; i++)
        {
            precursor += rkParams.lambda_i[i] * C[n][i] 
                / (1 + rkParams.lambda_i[i] * dt);
        }

        for(int g=0; g<G; g++)
        {
            S[index(n,g)] = phi[index(n,g)] * dx / (rkParams.v[g] * dt)
                    + rkParams.chi_d[g] * precursor * dx;
        }
    }
    return S;
}
   

//TODO: write description
Sparse formSigSMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape, 
        std::vector<double> adjoint)
{
    Sparse SigS = Sparse(mesh._G, mesh._G);

    for(int g=0; g < mesh._G; g++)
    {
        double diag = 0;
        for(int gp=0; gp < mesh._G; gp++)
        {
            if(gp != g)
            {
                double offdiag = 0;
                for(int n=0; n < mesh._N; n++)
                {
                    // calcualte diagonal element
                    XSdata* mat = mesh._material[n];
                    diag += mat->sigs[g][gp] * shape[index(n,g)] 
                        * adjoint[index(n,g)] * mesh._delta[n];
                    
                    // calculate off diagonal elements
                    offdiag += mat-> sigs[gp][g] * adjoint[index(n,g)]
                        * shape[index(n,gp)] * mesh._delta[n];
                }
                
                // set offdiagonal elements
                if(offdiag != 0)
                    SigS.setVal(g,gp,-offdiag);
            }
        }
        
        // set diagonal element
        if(diag != 0)
            SigS.setVal(g,g,diag);
                

    }
    return SigS;
}

//TODO: write description
Sparse formSigAMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape, 
        std::vector<double> adjoint)
{
    Sparse SigA = Sparse(mesh._G, mesh._G);
    
    for(int g=0; g < mesh._G; g++)
    {
        double val = 0;
        for(int n=0; n < mesh._N; n++)
        {
            // calcualte diagonal element
            XSdata* mat = mesh._material[n];
            val += mat->siga[g] * shape[index(n,g)] 
                * adjoint[index(n,g)] * mesh._delta[n];
        }

        // set diagonal element
        SigA.setVal(g,g,val);
    }
    return SigA;
}

// TODO: write description
Sparse formDhatMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape, 
        std::vector<double> adjoint)
{

    // re-create Dhat Matrix
    Sparse Dhat_SS = formDhatMatrix(mesh, index);
    
    // intialize condensed Dhat PKE Matrix
    Sparse Dhat = Sparse(mesh._G, mesh._G);

    // calculate collapsed diffusion coefficients
    for(int g=0; g < mesh._G; g++)
    {
        double val = 0;

        // sum over all mesh points
        for(int n=0; n < mesh._N; n++)
        {
            val += Dhat_SS(index(n,g),index(n,g)) 
                * shape[index(n,g)] * adjoint[index(n,g)]; 

            if( n != 0 )
                val += Dhat_SS(index(n,g), index(n-1,g)) * shape[index(n-1,g)]
                    * adjoint[index(n,g)];

            if( n != mesh._N-1 )
                val += Dhat_SS(index(n,g), index(n+1,g)) * shape[index(n+1,g)] 
                    * adjoint[index(n,g)];
        }

        Dhat.setVal(g,g,val);
    }

    return Dhat;
}

// TODO: write description
Sparse formFMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape,
       std::vector<double> adjoint, RKdata rkParams, double kcrit, double dt)
{
    Sparse F = Sparse(mesh._G, mesh._G);

    for(int g=0; g < mesh._G; g++)
    {
        for(int gp=0; gp < mesh._G; gp++)
        {
            double val = 0;
            for(int n=0; n < mesh._N; n++)
            {
                XSdata * mat = mesh._material[n];

                double prompt = (1 - rkParams.beta) * mat->chi[g];
                double delayed = 0;
                for(int i=0; i < rkParams.I; i++)
                    delayed += rkParams.beta_i[i] * rkParams.lambda_i[i] * dt
                        * rkParams.chi_d[g] / (1 + rkParams.lambda_i[i] * dt);

                val += (prompt + delayed) / kcrit * adjoint[index(n,g)]
                    * shape[index(n,gp)] * mat->nuSigf[gp] * mesh._delta[n];
            }

            // set fission matrix value
            if(val != 0)
                F.setVal(g,gp,val);
        }
    }
    return F;
}

//TODO: write description
std::vector<double> formSVectorPKE(Indexer index, std::vector<double> power, 
        RKdata rkParams, std::vector<std::vector<double> > C_tilde, double dt,
        std::vector<double> time_abs)
{
    std::vector<double> S (index._G, 0);

    for(int g=0; g < index._G; g++)
    {
        double val = power[g] * time_abs[g] / dt;
        for(int i=0; i < rkParams.I; i++)
        {
            val += rkParams.lambda_i[i] * C_tilde[g][i]
                / (1 + rkParams.lambda_i[i] * dt);
        }

        // assign value
        S[g] = val;
    }

    return S;
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

// sets the mesh and time steps for the transient
Transient::Transient()
{
    // set length
    n_pts = 0;
    set = false;
    tolerance = pow(10,-8);
}

// destructor for transient
Transient::~Transient() { }

// set interpolation times
void Transient::setInterpTimes(double * timeArray, int n_steps)
{
    timeVector.clear();
    n_pts = n_steps;
    for(int m=0; m < n_steps; m++)
        timeVector.push_back(timeArray[m]);
    return;
}
// set calculation times
void Transient::setCalcTimes(double * timeArray, int n_steps)
{
    timeSteps.clear();
    for(int m=0; m < n_steps; m++)
        timeSteps.push_back(timeArray[m]);
    return;
}

// set mesh
void Transient::setMeshVector(Mesh ** meshArray, int n_interp)
{
    if(n_interp != n_pts)
    {
        std::cout << "Error: number of mesh does not equal number of "
            << "interpolation times!" << endl;
        return;
    }
    meshVector.clear();
    for(int m=0; m < n_interp; m++)
        meshVector.push_back(*meshArray[m]);
    
    set = true;
    return;
}

// set tolerance
void Transient::setTolerance(double tol)
{
    tolerance = tol;
    return;
}
