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

rkSolution solvePKE(Transient trans, RKdata rkParams, int shape_step,
        int adj_weighting, int adj_step);

rkSolution solvePKESimple(Transient trans, RKdata rkParams, int shape_step,
        int adj_weighting, int adj_step);

rkSolution solveTransientFT(Transient trans, RKdata rkParams, int omegaMode);

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

Sparse formFhatMatrixFT(Sparse F, Mesh mesh, RKdata rkParams, double dt, 
        double kcrit, std::vector<double> omega, Indexer index);

std::vector<double> formSVectorFT(Mesh mesh, RKdata rkParams, 
        std::vector<double> source, std::vector<double> flux, 
        std::vector<std::vector<double> > C, double dt, double kcrit, 
        std::vector<double> omega, Indexer index);
#endif
