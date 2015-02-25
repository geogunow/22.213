#include"utils.h"

std::vector<double> zeros(int N)
{
    // initialize vector
    std::vector<double> * x = new std::vector<double>;

    // fill with zeros
    for(int i=0; i<N; i++)
        x->push_back(0);
    
    return *x;
}

