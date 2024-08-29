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

/*!
get the norm of the force vector
*/
double updater::getForceNorm()
    {
    vector<double> forceNorm(1);
    forceNorm[0] = 0.0;
    for (int ii = 0; ii < Ndof; ++ii)
        {
        forceNorm[0] += model->forces[ii].squared_length();
        }
    //define a lambda which is just addition
    sim->manipulateUpdaterData(forceNorm,
                         [](double x, double y)-> double
                                        {
                                        return x+y;
                                        });
    squaredTotalForceNorm = forceNorm[0];
    return sqrt(squaredTotalForceNorm);
    };

/*!
get the max norm of the force vector
*/
double updater::getMaxForce()
    {
    vector<double> maxNorm(1);
    maxNorm[0] = 0.0;

    for (int ii = 0; ii < Ndof; ++ii)
        {
        double currentNormSquared = model->forces[ii].squared_length();
        if (currentNormSquared > maxNorm[0])
            maxNorm[0] = currentNormSquared;
        }
    //define a lambda which is just the max operation
    sim->manipulateUpdaterData(maxNorm,
                         [](double x, double y)-> double
                                        {
                                        return std::max(x,y);
                                        });
    maximumForceNorm = sqrt(maxNorm[0]);
    return maximumForceNorm;
    };
