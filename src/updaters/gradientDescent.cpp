#include "gradientDescent.h"

void gradientDescent::performUpdate()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
    sim->computeForces();
    for (int ii = 0; ii < Ndof; ++ii)
        {
        displacements[ii] = deltaT*model->forces[ii];
        }
    sim->moveParticles(displacements);
    };

double gradientDescent::getForceNorm()
    {
    sim->computeForces();
    vector<double> forceNorm(1);
    forceNorm[0] = 0.0;
    for (int ii = 0; ii < Ndof; ++ii)
        {
        double currentNormSquared = model->forces[ii].squared_length();
        forceNorm[0] +=currentNormSquared;
printf("f %f\n",currentNormSquared);
        }
printf("f %f\n",forceNorm[0]);
    //define a lambda which is just addition
    sim->manipulateUpdaterData(forceNorm,
                         [](double x, double y)-> double 
                                        {
                                        return x+y;
                                        });
    return forceNorm[0];
    };
double gradientDescent::getMaxForce()
    {
    sim->computeForces();
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
    return sqrt(maxNorm[0]);
    };
