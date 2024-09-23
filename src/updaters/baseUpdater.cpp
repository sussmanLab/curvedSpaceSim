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
 get the norm of the force vector while excluding a subset of the particles, marked by particle index
 */
double updater::getForceNormWithExclusions(vector<int> exclusions)
    {
    vector<double> forceNorm(1);
    forceNorm[0] = 0.0;
    for (int ii = 0; ii < Ndof; ++ii)
        {
        if(find(exclusions.begin(), exclusions.end(), ii) != exclusions.end())
            {/* checks if exclusions contains ii -- we can cut down on computation time
                creating exclusions sorted, as indices are unique, then only ever
                considering the 'next' value of exclusions */
            continue;
            }
        forceNorm[0] += model->forces[ii].squared_length();
        }
    //define a lambda which is just addition
    sim->manipulateUpdaterData(forceNorm,
                         [](double x, double y)-> double
                                        {
                                        return x+y;
                                        });
    squaredTotalForceNorm = forceNorm[0];
    return squaredTotalForceNorm;
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


/*!
get the max norm of the force vector while excluding particle indices contained in exclusions
*/
double updater::getMaxForceWithExclusions(vector<int> exclusions)
    {
    vector<double> maxNorm(1);
    maxNorm[0] = 0.0;

    for (int ii = 0; ii < Ndof; ++ii)
        {
	if(find(exclusions.begin(), exclusions.end(), ii) != exclusions.end())
            {/* checks if exclusions contains ii -- we can cut down on computation time
                creating exclusions sorted, as indices are unique, then only ever
                considering the 'next' value of exclusions */
            continue;
            }
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
