#include"MonteCarlo.h"

/*
Function create a random number uniformly between 0 and 1
*/
double urand()
{
    return (double) rand() / RAND_MAX;
}


/*
Function to sample a unit vector assuming isotropic
*/
vec3D sample_flight_path()
{
    // initialize unit vector
    vec3D unit;

    // sample angles
    double phi = 2 * M_PI * urand();
    double mu = 2 * urand() - 1;

    // calculate unit vectors
    unit['x'] = sqrt(1 - pow(mu,2)) * cos(phi);
    unit['y'] = sqrt(1 - pow(mu,2)) * sin(phi);
    unit['z'] = mu;

    return unit;
}

/*
Creates initial neutron sites assuming an isotropic and uniform source
*/
std::stack<double> gen_initial_neutrons(ProblemDef P)
{
    srand(10);
    // initialize stack
    std::stack<double> initial;

    // loop over number of initial neutron histories
    for(int i=0; i < P.nhist; i++)
    {
        // add neutron to stack
        initial.push(0);
    }
    return initial;
}
/*
Creates initial precurosr sites assuming an isotropic and uniform distiribution
*/
std::stack<double> gen_initial_precursors(ProblemDef P)
{
    // initialize stack
    std::stack<double> precursors;

    // calculate number of precursors
    int n_precursors = P.beta * P.nhist * P.velocity * P.siga / P.lambda;

    // create precursors
    for(int i=0; i < n_precursors; i++)
    {
        // add neutron to stack
        precursors.push(0);
    }
    return precursors;
}

/*
Trace a neutron through a geometry, forming fission sites if encountered
*/
void trace_neutron(double &neutron, std::stack<double> &prompt, 
        std::stack<double> &delayed, int &stopped, ProblemDef P)
{
    // initialize vectors
    vec3D dist;
    std::vector<char> ref_surf;

    // set boolean tracker to alive
    bool alive = true;

    // follow neutron while alive
    while(alive)
    {
        // sample distance to collision
        double s = -log(urand()) / P.sigt;
        
        // determine if max distance in time crossed
        if(s + neutron > 1.0 / P.siga)
        {
            // add neutron to stopped bank
            stopped += 1;
            return;
        }

        // add distance traveled to the time distance
        neutron += s;

        // determine if neutron is absorbed in collision
        if(urand() < P.siga / P.sigt)
        {
            // kill current neutron
            alive = false;

            // determine if fission event
            if(urand() < P.sigf / P.siga)
            {
                // sample the number of neutrons created (2 or 3)
                int born;
                if(urand() < P.nu-2)
                    born = 3;
                else
                    born = 2;

                // create new neutrons (add to stack)
                for(int i=0; i < born; i++)
                {
                    // create neutron
                    double new_neutron = neutron;
                    
                    // decide whether to emit neutron prompt or delayed
                    if(urand() < P.beta)
                        delayed.push(new_neutron);
                    else
                        prompt.push(new_neutron);
                }
            }
        }
    }
    return;
}


/*
Follows neutrons in starting stack until completion of time step
*/
void follow_neutrons(int &starting, std::stack<double> &precursors, 
        ProblemDef P)
{
    // initialize prompt, delayed, and stopped stacks
    std::stack<double> prompt, delayed;
    int stopped = 0;
    
    // hand precursors to delayed stack
    delayed.swap(precursors);
        
    // initialize tracked neutron
    double neutron;

    // cycle through all stacks until all are empty
    while(!prompt.empty() or !delayed.empty() or starting != 0)
    {
        // pick a neutron
        if(starting != 0)
        {
            starting -= 1;
            neutron = 0;
        }
        else if(!prompt.empty())
        {
            neutron = prompt.top();
            prompt.pop();
        }
        else
        {
            // pick out a delayed neutron
            neutron = delayed.top();
            delayed.pop();

            // decide whether to use it or add to the stack
            double dt = (1.0/P.siga - neutron) / P.velocity;
            double prob_survival = exp(-P.lambda * dt);
            if(urand() < prob_survival)
            {
                neutron = 0;
                precursors.push(neutron);
                continue;
            }
            else
            {
                neutron += urand() * dt;
            }
        }

        // trace neturon through problem until death or time step completed
        trace_neutron(neutron, prompt, delayed, stopped, P);
    }
    starting = stopped;
}

        

