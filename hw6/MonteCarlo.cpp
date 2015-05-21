#include"MonteCarlo.h"

/*
Function create a random number uniformly between 0 and 1
*/
double urand()
{
    return (double) rand() / RAND_MAX;
}


/*
Trace a neutron through a geometry, forming fission sites if encountered
*/
void trace_neutron(double &neutron, std::stack<double> &prompt, 
        int &precursors, int &stopped, ProblemDef P)
{
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
                    {
                        // determine whether neutron decays in this timestep
                        double dt = (1.0/P.siga - neutron) / P.velocity;
                        double prob_survival = exp(-P.lambda * dt);
                        if(urand() < prob_survival)
                            precursors += 1;
                        else
                        {
                            new_neutron += urand() * dt;
                            prompt.push(new_neutron);
                        }

                    }
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
void follow_neutrons(int &starting, int &precursors, ProblemDef P)
{
    // initialize prompt and stopped stacks
    std::stack<double> prompt;
    int stopped = 0;
    
    // hand precursors to delayed stack
    int delayed = precursors;
    precursors = 0;
        
    // initialize tracked neutron
    double neutron;

    // cycle through all stacks until all are empty
    while(!prompt.empty() or delayed != 0 or starting != 0)
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
            delayed -= 1;

            // determine if delayed neutron will decay in timestep
            double dt = 1.0 / (P.siga * P.velocity);
            double prob_survival = exp(-P.lambda * dt);
            if(urand() < prob_survival)
            {
                precursors += 1;
                continue;
            }
            else
                neutron = urand() * dt;
        }

        // trace neturon through problem until death or time step completed
        trace_neutron(neutron, prompt, precursors, stopped, P);
    }
    starting = stopped;
}

        

