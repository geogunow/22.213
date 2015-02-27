#include"diffSolver.h"

int main(){
    // initialize fuel cross sections
    XSdata fuel1 = XSdata();
    int G = 2;
    double D[] = {1.3,0.5};
    double siga[] = {0.0098, 0.114};
    double nuSigf[] = {0.006, 0.195};
    double chi[] = {1.0, 0.0};
    double sigs_vals[] = {0,0.022,0,0};
    
    fuel1.setD(D,2);
    fuel1.setSiga(siga,2);
    fuel1.setNuSigf(nuSigf,2);
    fuel1.setChi(chi,2);
    for(int g=0; g<fuel1.G; g++)
    {
        std::vector<double> temp;
        for(int gp=0; gp<fuel1.G; gp++)
        {
            temp.push_back(sigs_vals[fuel1.G*g + gp]);
        }
        fuel1.sigs.push_back(temp);
    }

    // initialize water cross sections
    XSdata water = XSdata();
    double D2[] = {1.3,0.5};
    double siga2[] = {0.0002, 0.010};
    double nuSigf2[] = {0.0, 0.0};
    double chi2[] = {1.0, 0.0};
    double sigs_vals2[] = {0,0.032,0,0};
    
    water.setD(D2,2);
    water.setSiga(siga2,2);
    water.setNuSigf(nuSigf2,2);
    water.setChi(chi2,2);
    for(int g=0; g<fuel1.G; g++)
    {
        std::vector<double> temp;
        for(int gp=0; gp<water.G; gp++)
        {
            temp.push_back(sigs_vals2[water.G*g + gp]);
        }
        water.sigs.push_back(temp);
    }

    // form mesh
    int nregions = 5;
    XSdata* mats[] = {&water, &fuel1, &fuel1, &fuel1, &water};
    double widths[] = {20.0, 15.0, 20.0, 15.0, 20.0};
    int npts[] = {20, 20, 30, 20, 20};
    
    Mesh mesh;
    mesh.setMesh(npts, nregions);
    mesh.setWidths(widths, nregions);
    mesh.setMaterials(mats, nregions);
    mesh.setBoundaryConditions(0, 0);

    // solve eigenvalue problem
    double outer_tol = pow(10,-6);
    double inner_tol = pow(10,-6);
    int maxiters = 10000;
    int inner_solver = 2; // gauss seidel (SOR)
    eigenSolution result = solveDiffusion(mesh, outer_tol, inner_tol,
            maxiters, inner_solver);
   
    std::vector<double> F = result.getFlux(1);
    /*
    for(int i=0; i<F.size(); i++)
        std::cout << F[i] << endl;
    return 0;
    */
}

