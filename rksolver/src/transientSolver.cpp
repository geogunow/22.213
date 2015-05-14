#include"transientSolver.h"

//TODO: Description
rkSolution solveTransient(Transient trans, RKdata rkParams)
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

// sets the mesh and time steps for the transient
Transient::Transient()
{
    // set length
    n_pts = 0;
    set = false;
    tolerance = pow(10,-6);
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
