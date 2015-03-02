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
    //int npts[] = {20, 20, 30, 20, 20};
    int npts[] = {1, 1, 1, 1, 1};
    
    Mesh mesh;
    mesh.setMesh(npts, nregions);
    mesh.setWidths(widths, nregions);
    mesh.setMaterials(mats, nregions);
    mesh.setBoundaryConditions(0, 0);

    // solve eigenvalue problem
    double outer_tol = pow(10,-6);
    double inner_tol = pow(10,-6);
    int maxiters = 10000;
    int inner_solver = 2; // optimal SOR
    eigenSolution result = solveCritical(mesh, outer_tol, inner_tol,
            maxiters, inner_solver);
   
    std::vector<double> F = result.getFlux(1);
    std::cout << "keff = " << result.keff << endl;
    
    // transient solver
    std::cout << "Calculating transient...." << endl;

    std::vector<Mesh> meshVector;
    std::vector<double> timeVector;
    std::vector<double> timeSteps;

    double tt[] = {0, 1, 5};
    for(int i=0; i<3; i++)
    {
        timeVector.push_back(tt[i]);
        meshVector.push_back(mesh);
    }
    double ts[] = {0,1,2,3,4};
    for(int i=0; i<5; i++)
        timeSteps.push_back(ts[i]);

    int I = 3;
    double chi_d[] = {0,1};
    double v[] = {1000, 0.001};
    double beta[] = {0.0030, 0.0030, 0.0005};
    double lambda[] = {2.4, 1.2, 0.7};

    
    RKdata rkParams;
    rkParams.setChiD(chi_d, G);
    rkParams.setV(v, G);
    rkParams.setBetas(beta, I);
    rkParams.setLambdas(lambda, I);
    
    
    solveTransient( meshVector, timeVector, timeSteps, rkParams);
    std::cout << "Success!" << endl;
    return 0;
}

