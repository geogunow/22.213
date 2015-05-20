#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stack>
#include<array>
#include<vector>
#include<math.h>

#define endl std::endl

// custom types and data structures
class vec3D
{
    private:
    std::array<double,3> data;
    
    public:
    double& operator[](char var)
    {
        return data[(int) var - 120];
    }
    const double& operator[](char var) const
    {
        return data[(int) var - 120];
    }
};

struct ProblemDef
{
    double width;
    double nu;
    double sigs;
    double sigc;
    double sigf;
    double siga;
    double sigt;
    double beta;
    double lambda;
    double velocity;
    int ngen;
    int nhist;
};

struct Neutron
{
    double time_dist;
};

// functions
double urand();
vec3D sample_flight_path();
std::stack<double> gen_initial_neutrons(ProblemDef P);
std::stack<double> gen_initial_precursors(ProblemDef P);
void trace_neutron(double &neutron, std::stack<double> &prompt, 
        std::stack<double> &delayed, std::stack<double> &stopped, 
        ProblemDef P);
void follow_neutrons(int &starting, std::stack<double> &precursors, 
        ProblemDef P);
