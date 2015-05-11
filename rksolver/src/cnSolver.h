#ifndef __CNSOLVER__
#define __CNSOLVER__
#include"transientSolver.h"

rkSolution solveTransientCN(Transient trans, RKdata rkParams);
Sparse formFhatMatrixCN(Mesh mesh, RKdata rkParams, double dt, double kcrit, 
        Indexer index);
std::vector<double> formSVectorCN(Mesh mesh, RKdata rkParams, 
        std::vector<double> phi, std::vector<double> phi_NB,
        std::vector<std::vector<double> > C, double dt, double kcrit, 
        Indexer index);

#endif
