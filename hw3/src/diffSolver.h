#ifndef __DIFFSOLVER__
#define __DIFFSOLVER__
#include"utils.h"
#include"Sparse.h"
#include"Mesh.h"
#include<cmath>


class eigenSolution{
    public:
        std::vector<double> flux;
        std::vector< std::vector<double> > gFlux;
        std::vector<double> power;
        std::vector< std::vector<double> > gPower;
        double keff;
        int outer_iters;
        int inner_iters;
        
        eigenSolution();
        virtual ~eigenSolution();
        void indexArrays(Indexer index);
        std::vector<double> getFlux(int g);
        std::vector<double> getPower(int g);
};


eigenSolution solveCritical(Mesh mesh, double outer_tol, double inner_tol,
        int maxiters, int inner_solver);

void solveTransient(
        std::vector<Mesh> meshVector,
        std::vector<double> timeVector,
        std::vector<double> timeSteps,
        RKdata rkParams);

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
