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
std::stack<Neutron> gen_initial_neutrons(ProblemDef P)
{
    // initialize stack
    std::stack<Neutron> initial;

    // loop over number of initial neutron histories
    for(int i=0; i < P.nhist; i++)
    {
        // initialize neutron
        Neutron neutron;

        // initialize position
        char vars[] = {'x', 'y', 'z'};
        for(int j=0; j < 3; j++)
        {
            char var = vars[j];
            neutron.loc[var] = urand() * P.width;
        }

        // initialize direction
        neutron.unit = sample_flight_path();

        // add neutron to stack
        initial.push(neutron);
    }
    return initial;
}


/*
Trace a neutron through a geometry, forming fission sites if encountered
*/
void trace_neutron(Neutron &neutron, std::stack<Neutron> &prompt, 
        std::stack<Neutron> &delayed, std::stack<Neutron> &stopped, 
        ProblemDef P)
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
        
        // complete track length s through reflections
        bool tracking = true;
        while(tracking)
        {
            // calculate distance to intersection
            char vars[] = {'x', 'y', 'z'};
            for(int i=0; i < 3; i++)
            {
                char var = vars[i];
                double pos = (P.width - neutron.loc[var]) / neutron.unit[var];
                double neg = -neutron.loc[var] / neutron.unit[var];
                dist[var] = std::max(pos, neg);
            }

            // calculate min distance
            double min_dist = dist['x'];
            for(int i=0; i < 3; i++)
            {
                char var = vars[i];
                if(dist[var] < min_dist)
                {
                    ref_surf.clear();
                    ref_surf.push_back(var);
                    min_dist = dist[var];
                }
                else if(dist[var] == min_dist)
                    ref_surf.push_back(var);
            }

            // shorten distance
            s -= min_dist;

            // determine if max distance in time crossed
            if(min_dist + neutron.time_dist > 1.0 / P.siga)
            {
                // move neutron
                for(int i=0; i < 3; i++)
                {
                    char var = vars[i];
                    neutron.loc[var] += neutron.unit[var] * 
                        (1.0 / P.siga - neutron.time_dist);
                }
                
                // add neutron to stopped bank
                stopped.push(neutron);

                return;
            }

            // move neutron
            for(int i=0; i < 3; i++)
            {
                char var = vars[i];
                neutron.loc[var] += neutron.unit[var] * min_dist;
            }
          

            // add distance traveled to the time distance
            neutron.time_dist += min_dist;

            // reflect off surfaces if applicable
            for(int i=0; i < 3; i++)
            {
                char var = vars[i];
                neutron.unit[var] *= -1;
            }

            // check if original track length s is completed
            if(s <= 0)
                tracking = false;
        }

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
                    Neutron new_neutron;
                    new_neutron.loc = neutron.loc;
                    new_neutron.unit = sample_flight_path();
                    new_neutron.time_dist = neutron.time_dist;
                    
                    // decide whether to emit neutron prompt or delayed
                    if(urand() < P.beta)
                        delayed.push(new_neutron);
                    else
                        prompt.push(new_neutron);
                }
            }
        }
        
        // scatter neutron
        else
            neutron.unit = sample_flight_path();
    }
    return;
}


/*
Follows neutrons in starting stack unitl completion of time step
*/
void follow_neutrons(std::stack<Neutron> &starting, 
        std::stack<Neutron> &precursors, ProblemDef P)
{
    // initialize prompt, delayed, and stopped stacks
    std::stack<Neutron> prompt, delayed, stopped;
    
    // hand precursors to delayed stack
    delayed.swap(precursors);

    // cycle through all stacks until all are empty
    while(!prompt.empty() or !delayed.empty() or !starting.empty())
    {
        // pick a neutron
        Neutron neutron;
        if(!starting.empty())
        {
            neutron = starting.top();
            starting.pop();
            neutron.time_dist = 0;
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
            double dt = (1.0/P.siga - neutron.time_dist) / P.velocity;
            double prob_survival = exp(-P.lambda * dt);
            if(urand() < prob_survival)
            {
                neutron.time_dist = 0;
                precursors.push(neutron);
                continue;
            }
        }

        // trace neturon through problem until death or time step completed
        trace_neutron(neutron, prompt, delayed, stopped, P);
    }
    starting = stopped;
}

        

