#ifndef __DIFFSOLVER__
#define __DIFFSOLVER__
#include"utils.h"
#include"Sparse.h"
#include"Mesh.h"
#include"Solutions.h"
#include<cmath>


class Transient{
    public:
        int n_pts;
        bool set;
        std::vector<Mesh> meshVector;
        std::vector<double> timeVector;
        std::vector<double> timeSteps;

        Transient();
        virtual ~Transient();
        void setInterpTimes(double * timeArray, int n_steps);
        void setCalcTimes(double * timeArray, int n_steps);
        void setMeshVector(Mesh ** meshArray, int n_interp);
};


eigenSolution solveCritical(Mesh mesh, double outer_tol, double inner_tol,
        int maxiters, int inner_solver);

rkSolution solveTransient(Transient trans, RKdata rkParams);

Sparse formSigAMatrix(Mesh mesh, Indexer index);
Sparse formSigSMatrix(Mesh mesh, Indexer index);
Sparse formDhatMatrix(Mesh mesh, Indexer index);
Sparse formFMatrix(Mesh mesh, Indexer index);
Sparse formFhatMatrix(Mesh mesh, RKdata rkParams, double dt, double kcrit, 
        Indexer index);
std::vector<double> formSVector(Mesh mesh, RKdata rkParams, 
        std::vector<double> phi, std::vector<std::vector<double> > C, 
        double dt, double kcrit, Indexer index);
void formSteadyStateMatrixProblem(Mesh mesh, Sparse &A, Sparse &F, 
        Indexer index);

eigenSolution eigen_solver(Sparse A, Sparse F, int N, int G, double tol, 
        double inner_tol, int maxiters, int inner_solver);

#endif
