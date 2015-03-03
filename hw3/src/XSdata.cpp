#include"XSdata.h"
/*
   Constructor for mesh
   */
XSdata::XSdata()
{
    G = 0;
}

/*
   Destructor for mesh
   */
XSdata::~XSdata() { }


/* 
   functions for setting material properties
   */
void XSdata::setD(double * xs, int ng)
{
    if(G == 0){
        G = ng;
        if(ng == 0)
        {
            std::cout << "Error: groups set to 0" << endl;
            return;
        }
    }
    else if(ng != G)
    {
        std::cout << "Error group mismatch" << endl;
        return;
    }

    // fill vector with data
    for(int g=0; g<G; g++)
        D.push_back(xs[g]);
    return;
}
void XSdata::setSiga(double * xs, int ng)
{
    if(G == 0){
        G = ng;
        if(ng == 0)
        {
            std::cout << "Error: groups set to 0" << endl;
            return;
        }
    }
    else if(ng != G)
    {
        std::cout << "Error group mismatch" << endl;
        return;
    }

    // fill vector with data
    for(int g=0; g<G; g++)
        siga.push_back(xs[g]);
    return;
}
void XSdata::setNuSigf(double * xs, int ng)
{
    if(G == 0){
        G = ng;
        if(ng == 0)
        {
            std::cout << "Error: groups set to 0" << endl;
            return;
        }
    }
    else if(ng != G)
    {
        std::cout << "Error group mismatch" << endl;
        return;
    }

    // fill vector with data
    for(int g=0; g<G; g++)
        nuSigf.push_back(xs[g]);
    return;
}
void XSdata::setChi(double * xs, int ng)
{
    if(G == 0){
        G = ng;
        if(ng == 0)
        {
            std::cout << "Error: groups set to 0" << endl;
            return;
        }
    }
    else if(ng != G)
    {
        std::cout << "Error: group mismatch" << endl;
        return;
    }

    // fill vector with data
    for(int g=0; g<G; g++)
        chi.push_back(xs[g]);
    return;
}
void XSdata::setSigs(double ** xs, int ng1, int ng2)
{

    if(ng1 != ng2)
    {
        std::cout << "Error: scattering matrix group mismatch" << endl;
        return;
    }

    if(G == 0){
        G = ng1;
        if(ng1 == 0)
        {
            std::cout << "Error: groups set to 0" << endl;
            return;
        }
    }
    else if(ng1 != G)
    {
        std::cout << "Error: group mismatch" << endl;
        return;
    }

    // fill vector with data
    for(int g=0; g<G; g++)
    {
        std::vector<double> temp;
        for(int gp=0; gp<G; gp++)
            temp.push_back(xs[g][gp]);
        sigs.push_back(temp);
    }
    return;
}

void XSdata::printXS()
{
    for(int i=0; i<G; i++)
    {
        for(int j=0; j<G; j++)
        {
            if(j == 0)
                std::cout << "[";
            std::cout << sigs[i][j];
            if(j == G-1)
                std::cout << "]\n";
            else
                std::cout << ", ";
        }
    }
    return;
}

/*
   Constructor for reactor kinetics parameters
   */
RKdata::RKdata()
{
    I = 0;
    G = 0;
    beta = 0;
}

/*
   Destructor for reactor kinetics parameters
   */
RKdata::~RKdata() { }


/* 
   functions for setting reactor kinetics parameters
   */
void RKdata::setV(double * vals, int len)
{
    // error checking
    if( G == 0 )
    {
        G = len;
        if( len == 0 )
        {
            std::cout << "Error: groups set to 0" << endl;
            return;
        }
    }
    else if( len != G )
    {
        std::cout << "Error: group mismatch" << endl;
        return;
    }
    if( v.size() != 0 )
    {
        std::cout << "Error: velocities already set" << endl;
        return;
   }

    // set velocities
    for(int g=0; g < G; g++)
        v.push_back( vals[g] );

    return;
}
void RKdata::setChiD(double * vals, int len)
{
    // error checking
    if( G == 0 )
    {
        G = len;
        if( len == 0 )
        {
            std::cout << "Error: groups set to 0" << endl;
            return;
        }
    }
    else if( len != G )
    {
        std::cout << "Error: group mismatch" << endl;
        return;
    }
    if( chi_d.size() != 0 )
    {
        std::cout << "Error: delayed fission spectrum already set" << endl;
        return;
   }

    // set delayed fission spectrum
    for(int g=0; g < G; g++)
        chi_d.push_back( vals[g] );

    return;
}
void RKdata::setBetas(double * vals, int len)
{
    // error checking
    if( I == 0 )
    {
        I = len;
        if( len == 0 )
        {
            std::cout << "Error: precursor groups set to 0" << endl;
            return;
        }
    }
    else if( len != I )
    {
        std::cout << "Error: precursor group mismatch" << endl;
        return;
    }
    if( beta_i.size() != 0 )
    {
        std::cout << "Error: delayed neutron fractions already set" << endl;
        return;
    }

    // set delayed neutron fractions
    for(int i=0; i < I; i++)
    {
        beta += vals[i];
        beta_i.push_back( vals[i] );
    }

    return;
}
void RKdata::setLambdas(double * vals, int len)
{
    // error checking
    if( I == 0 )
    {
        I = len;
        if( len == 0 )
        {
            std::cout << "Error: precursor groups set to 0" << endl;
            return;
        }
    }
    else if( len != I )
    {
        std::cout << "Error: precursor group mismatch" << endl;
        return;
    }
    if( lambda_i.size() != 0 )
    {
        std::cout << "Error: delayed neutron decay already set" << endl;
        return;
   }

    // set delayed neutron decay constants
    for(int i=0; i < I; i++)
        lambda_i.push_back( vals[i] );

    return;
}



