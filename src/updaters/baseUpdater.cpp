#include "baseUpdater.h"
/*! \file baseUpdater.cpp" */

void updater::Update(int timestep)
    {
    iterations = timestep;
    if (maxIterations >0 && maxIterations < iterations)
        return;
    if(Period <= 0 || (Period >0 && (timestep+Phase) % Period == 0))
        performUpdate();
    };

void updater::performUpdate()
    {
    cout << "in the base updater... that's odd..." << endl;
    sim->computeForces();
    };
