#include"MonteCarlo.h"

int main()
{
    // define problem
    ProblemDef P;
    P.width = 400;
    P.nu = 2.5;
    P.sigs = 0.27;
    P.sigc = 0.02;
    P.sigf = P.sigc / (P.nu - 1);
    P.siga = P.sigf + P.sigc;
    P.sigt = P.sigs + P.siga;
    P.beta = 0; // 0.01
    P.lambda = log(2.0) / 10.0;
    P.velocity = 2200 * sqrt(0.1/0.0253) * 100; //0.1 eV neutrons assumed
    P.ngen = 1500;
    P.nhist = pow(10,4);

    // define transient
    int transient_step = 100;
    double transient_magnitude = 1.00; //1.002;

    // initialize vectors to store information
    std::vector<int> n_pop;
    std::vector<vec3D> n_locs;
   
    // generate initial neutrons
    std::stack<Neutron> starting = gen_initial_neutrons(P);
    std::stack<Neutron> precursors = gen_initial_precursors(P);

    // record initial neutron locations and population
    n_pop.push_back(P.nhist);

    // follow neutrons through generations
    for(int t=0; t < P.ngen; t++)
    {
        // display progress
        if( (t+1) % 100 == 0 )
            std::cout << "Step " << t+1 << "/" << P.ngen << endl;

        // transient definition
        if(t == transient_step-1)
            P.nu *= transient_magnitude;

        // follow neutrons to the end of the time step
        follow_neutrons(starting, precursors, P);

        // record neutron locations and populations at each timestep
        n_pop.push_back(starting.size());
    }

    // print results to output file
    std::ofstream out;
    out.open("out.txt");
    for(int t=0; t < n_pop.size(); t++)
        out << n_pop[t] << endl;
    out.close();
    
    return 0;
}

