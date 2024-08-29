#include "velocityVerletNVE.h"

void velocityVerletNVE::performUpdate()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
    //steps 1 and 2: update next displacement set, and do first half of the velocity update
    velocityVerletFirstHalfStep();
    //step 3a: move particles and compute the new forces
    //step 3b: with the new forces, compute the second half of the velocity update
    velocityVerletSecondHalfStep();
    };

void velocityVerletNVE::velocityVerletFirstHalfStep()
    {
    for (int ii = 0; ii < Ndof; ++ii)
        {
        displacements[ii] = deltaT*model->velocities[ii] + 0.5*deltaT*deltaT*model->forces[ii];
        model->velocities[ii] += 0.5*deltaT*model->forces[ii];
        }
    };

void velocityVerletNVE::velocityVerletSecondHalfStep()
    {
    sim->moveParticles(displacements);
    sim->computeForces();
    for (int ii = 0; ii < Ndof; ++ii)
        model->velocities[ii] += 0.5*deltaT*model->forces[ii];
    };
