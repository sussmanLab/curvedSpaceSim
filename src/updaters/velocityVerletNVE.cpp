#include "velocityVerletNVE.h"

void velocityVerletNVE::performUpdate()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
    //steps 1 and 2: update next displacement set, and do first half of the velocity update
    for (int ii = 0; ii < Ndof; ++ii)
        {
        displacements[ii] = deltaT*model->velocities[ii] + 0.5*deltaT*deltaT*model->forces[ii];
        model->velocities[ii] += 0.5*deltaT*model->forces[ii];
        }
    //step 3a: move particles and compute the new forces
    sim->moveParticles(displacements);
    sim->computeForces();
    //setp 3b: with the new forces, compute the second half of the velocity update
    for (int ii = 0; ii < Ndof; ++ii)
        model->velocities[ii] += 0.5*deltaT*model->forces[ii];
    };
