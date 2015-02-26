#include"diffSolver.h"

int main(){

    /*
    Sparse A = Sparse(3,3);
    A.setVal(0,0,5);
    A.setVal(0,1,10); 
    A.setVal(1,1,4); 

    A.setVal(2,2,3); 
    std::vector<std::vector<double> > Adense;
    Adense = A.dense();

    A.display();
    
    Indexer index = Indexer(2,2);
    
    Sparse T = A.transpose();

    std::cout << endl;

    T.display();
    
    // test matrix multiply
    std::vector<double> xx (3,1);
    xx[1] = 10;
    std::vector<double> b = T*xx;
    for(int i=0; i<3; i++)
        std::cout << b[i] << endl;
        */
    
    // initialize cross sections
    XSdata fuel1;
    int G = 2;
    double D[] = {1.3,0.5};
    double siga[] = {0.0098, 0.114};
    double nuSigf[] = {0.006, 0.195};
    double chi[] = {1.0, 0.0};
    double sigs_vals[] = {0,0.022,0,0};
    
    fuel1.D = std::vector<double> (D, D + sizeof(D) / sizeof(double));
    fuel1.siga = std::vector<double> (siga, siga + sizeof(siga) / sizeof(double));
    fuel1.nuSigf = std::vector<double> 
        (nuSigf, nuSigf + sizeof(nuSigf) / sizeof(double));
    fuel1.chi = std::vector<double> (chi, chi + sizeof(chi) / sizeof(double));
    for(int g=0; g<G; g++)
    {
        std::vector<double> temp;
        for(int gp=0; gp<G; gp++)
        {
            temp.push_back(sigs_vals[G*g + gp]);
        }
        fuel1.sigs.push_back(temp);
    }

    // form mesh
    XSdata* mats[] = {&fuel1, &fuel1, &fuel1};
    double w[] = {1.5, 2.0, 1.5};
    int npt[] = {2, 3, 2};
    
    std::vector<XSdata*> materials (mats, mats + sizeof(mats) / sizeof(XSdata*));
    std::vector<double> widths (w, w + sizeof(w) / sizeof(double));
    std::vector<int> npts (npt, npt + sizeof(npt) / sizeof(int));

    Mesh mesh = fillGeometry(materials, widths, npts);
    
    // form matrices
    int N = mesh.delta.size();
    BoundaryConditions BC;
    BC.left = 0;
    BC.right = 0;
    Sparse A = Sparse(N*G, N*G);
    Sparse F = Sparse(N*G, N*G);
    Indexer index = Indexer(N,G);
    
    formMatrixProblem(mesh, BC, A, F, index);

    //F.display();

    double outer_tol = pow(10,-6);
    double inner_tol = pow(10,-6);
    int maxiters = 1000;
    int inner_solver = 1; // gauss seidel
    eigenSolution result = eigen_solver(A, F, N, G, outer_tol, inner_tol,
            maxiters, inner_solver);
    std::cout << "keff = " << result.keff << endl;
    
    return 0;
}

