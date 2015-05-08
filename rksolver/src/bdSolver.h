#ifndef __BDSOLVER__
#define __BDSOLVER__
#include"transientSolver.h"

rkSolution solveTransientBD(Transient trans, RKdata rkParams, int BD);
Sparse formFhatMatrixBD(Mesh mesh, RKdata rkParams, double dt, double kcrit, 
        Indexer index, double alpha, double omega);
std::vector<double> formSVectorBD(Mesh mesh, RKdata rkParams, 
        std::vector<double> F_sum, std::vector<std::vector<double> > C_sum,
        double dt, double kcrit, Indexer index, double alpha, double omega);

#endif
