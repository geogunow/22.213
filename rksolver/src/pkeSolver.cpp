#include"pkeSolver.h"

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
