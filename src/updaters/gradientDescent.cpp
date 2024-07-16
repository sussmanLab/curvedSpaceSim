#include "gradientDescent.h"

void gradientDescent::performUpdate()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
    sim->computeForces();
    for (int ii = 0; ii < Ndof; ++ii)
        {
        displacements[ii] = deltaT*model->forces[ii];
        //printf("p %i d (%f %f %f)  %f\n",ii,displacements[ii][0],displacements[ii][1],displacements[ii][2], displacements[ii].squared_length());
        }
    sim->moveParticles(displacements);
    };

double gradientDescent::getForceNorm()
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
    return squaredTotalForceNorm;
    };

double gradientDescent::getMaxForce()
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
