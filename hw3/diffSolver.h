#include"utils.h"
#include"Sparse.h"
#include<cmath>

typedef struct{
    std::vector<double> D;
    std::vector<double> siga;
    std::vector<std::vector<double> > sigs;
    std::vector<double> nuSigf;
    std::vector<double> chi;
} XSdata;

typedef struct{
    std::vector<double> delta;
    std::vector<XSdata*> material;
} Mesh;

typedef struct{
    int left;
    int right;
} BoundaryConditions;

typedef struct{
    std::vector<double> flux;
    std::vector<double> power;
    double keff;
    int outer_iters;
    int inner_iters;
} eigenSolution;



class Indexer{
    public:
        int _N;
        int _G;
        Indexer(int N, int G)
        {
            _N = N;
            _G = G;
        };
        virtual ~Indexer() {};
        int operator() (int n, int g){
            return _N*g + n;
        };
};

void formMatrixProblem(Mesh mesh, BoundaryConditions BC, 
        Sparse &A, Sparse &F, Indexer index);
Mesh fillGeometry(std::vector<XSdata*> materials, std::vector<double> widths, 
        std::vector<int> npts);

eigenSolution eigen_solver(Sparse A, Sparse F, int N, int G, double tol, 
        double inner_tol, int maxiters, int inner_solver);

std::vector<double> pointJacobi(Sparse A, std::vector<double> b, double tol, 
        int maxiters, int &sumiters);

std::vector<double> gaussSeidel(Sparse A, std::vector<double> b, double tol, 
        int maxiters, int &sumiters);
