#include"utils.h"

/*
   Function to compute the dot product between two vectors
   */
double dot(std::vector<double> x1, std::vector<double> x2)
{
    // check vector sizes
    if(x1.size() != x2.size())
        std::cout << "Error: vector size mismatch" << endl;

    // calculate dot product
    double s = 0;
    for(int i=0; i < x1.size(); i++)
        s += x1[i] * x2[i];

    return s;
}

/*
   Function to display a vector of doubles
   */
void printVector(std::vector<double> x)
{
    std::printf("[");
    for(int i=0; i<x.size(); i++)
    {
        std::printf("%f", x[i]);
        if( i != x.size()-1 )
            std::printf(",\n");
        else
            std::printf("]\n");
    }
    return;
}

