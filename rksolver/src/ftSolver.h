#ifndef __FTSOLVER__
#define __FTSOLVER__
#include"transientSolver.h"

rkSolution solveTransientFT(Transient trans, RKdata rkParams, int omegaMode);

Sparse formFhatMatrixFT(Sparse F, Mesh mesh, RKdata rkParams, double dt, 
        double kcrit, std::vector<double> omega, Indexer index);

std::vector<double> formSVectorFT(Mesh mesh, RKdata rkParams, 
        std::vector<double> source, std::vector<double> flux, 
        std::vector<std::vector<double> > C, double dt, double kcrit, 
        std::vector<double> omega, Indexer index);

#endif
