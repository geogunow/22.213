#ifndef __PKESOLVER__
#define __PKESOLVER__
#include"transientSolver.h"

rkSolution solvePKE(Transient trans, RKdata rkParams, int shape_step,
        int adj_weighting, int adj_step);

rkSolution solvePKESimple(Transient trans, RKdata rkParams, int shape_step,
        int adj_weighting, int adj_step);

Sparse formSigSMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape, 
        std::vector<double> adjoint);

Sparse formSigAMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape, 
        std::vector<double> adjoint);

Sparse formDhatMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape, 
        std::vector<double> adjoint);

Sparse formFMatrixPKE(Mesh mesh, Indexer index, std::vector<double> shape,
       std::vector<double> adjoint, RKdata rkParams, double kcrit, double dt);

std::vector<double> formSVectorPKE(Indexer index, std::vector<double> power, 
        RKdata rkParams, std::vector<std::vector<double> > C_tilde, double dt,
        std::vector<double> time_abs);


#endif
