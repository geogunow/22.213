#include"Solutions.h"

/*
   Eigenvalue solution constructor
   */
eigenSolution::eigenSolution()
{
    keff = 0;
    outer_iters = 0;
    inner_iters = 0;
}

//  Eigenvalue solution destructor
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
   Returns the flux profile for group g
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
        return power;
    }
    else
        return gPower[g-1];
} 

/*
   Reactor kinetics solution constructor
   */
rkSolution::rkSolution()
{
    unused = 0;
}

//  Reactor kinetics solution destructor
rkSolution::~rkSolution() {};

/*
   creates 3D arrays of flux and power by timestep and group
   */
void rkSolution::indexArrays(Indexer index, int T)
{
    if( gFlux.size() > 0 )
    {
        std::cout << "Warning: arrays already indexed" << endl;
        return;
    }

    // fill arrays
    for(int t=0; t < T; t++)
    {
        std::vector<std::vector<double> > tempFlux1;
        std::vector<std::vector<double> > tempPower1;
        for(int g=0; g<index._G; g++)
        {
            std::vector<double> tempFlux2;
            std::vector<double> tempPower2;
            for(int n=0; n<index._N; n++)
            {
                tempFlux2.push_back( flux[t][index(n,g)] );
                tempPower2.push_back( powerProfile[t][index(n,g)] );
            }
            tempFlux1.push_back(tempFlux2);
            tempPower1.push_back(tempPower2);
        }
        gFlux.push_back(tempFlux1);
        gPowerProfile.push_back(tempPower1);
    }
    return;
}

/*
   Returns the flux profile at time step t for group g
   */       
std::vector<double> rkSolution::getFlux(int t, int g)
{
    if(gFlux.size() == 0)
    {
        std::cout << "Error: arrays not indexed!" << endl;
        return flux[t];
    }
    else
        return gFlux[t][g-1];
}

/*
   Function to return a vector of the fission reaction rate for group g
   */
std::vector<double> rkSolution::getPowerProfile(int t, int g)
{
    if(gPowerProfile.size() == 0)
    {
        std::cout << "Error: arrays not indexed!" << endl;
        return powerProfile[t];
    }
    else
        return gPowerProfile[t][g-1];
} 

/*
   Function to get raw power profile vector at time t
   */
std::vector<double> rkSolution::getRawPowerProfile(int t)
{
    return powerProfile[t];
}


/*
   Returns a vector of the total integrated power at each time step
   */
std::vector<double> rkSolution::getPower()
{
    return power;
}
