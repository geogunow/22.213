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


eigenSolution solveDiffusion(Mesh mesh, double outer_tol, double inner_tol,
        int maxiters, int inner_solver);

void formMatrixProblem(Mesh mesh, Sparse &A, Sparse &F, Indexer index);

eigenSolution eigen_solver(Sparse A, Sparse F, int N, int G, double tol, 
        double inner_tol, int maxiters, int inner_solver);

#endif
