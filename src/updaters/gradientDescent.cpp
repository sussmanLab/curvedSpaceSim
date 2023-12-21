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

double gradientDescent::getMaxForce()
    {
    double maximumNorm = 0;

    for (int ii = 0; ii < Ndof; ++ii)
        {
        double currentNormSquared = model->forces[ii].squared_length();
        if (currentNormSquared > maximumNorm)
            maximumNorm = currentNormSquared;
        }
    return sqrt(maximumNorm);
    };
