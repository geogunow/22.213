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
