#ifndef __DIFFSOLVER__
#define __DIFFSOLVER__
#include"utils.h"
#include"Sparse.h"
#include"Mesh.h"
#include"Solutions.h"
#include<cmath>

eigenSolution solveCritical(Mesh mesh, double outer_tol, double inner_tol,
        int maxiters, int inner_solver);

Sparse formSigAMatrix(Mesh mesh, Indexer index);
Sparse formSigSMatrix(Mesh mesh, Indexer index);
Sparse formDhatMatrix(Mesh mesh, Indexer index);
Sparse formFMatrix(Mesh mesh, Indexer index);
void formSteadyStateMatrixProblem(Mesh mesh, Sparse &A, Sparse &F, 
        Indexer index);

eigenSolution eigen_solver(Sparse A, Sparse F, int N, int G, double tol, 
        double inner_tol, int maxiters, int inner_solver);

#endif
