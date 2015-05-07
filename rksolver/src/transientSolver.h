#ifndef __TRANSSOLVER__
#define __TRANSSOLVER__
#include"diffSolver.h"

class Transient{
    public:
        int n_pts;
        bool set;
        double tolerance;
        std::vector<Mesh> meshVector;
        std::vector<double> timeVector;
        std::vector<double> timeSteps;

        Transient();
        virtual ~Transient();
        void setInterpTimes(double * timeArray, int n_steps);
        void setCalcTimes(double * timeArray, int n_steps);
        void setMeshVector(Mesh ** meshArray, int n_interp);
        void setTolerance(double tol);
};


rkSolution solveTransient(Transient trans, RKdata rkParams);
Sparse formFhatMatrix(Mesh mesh, RKdata rkParams, double dt, double kcrit, 
        Indexer index);
std::vector<double> formSVector(Mesh mesh, RKdata rkParams, 
        std::vector<double> phi, std::vector<std::vector<double> > C, 
        double dt, double kcrit, Indexer index);

#endif
