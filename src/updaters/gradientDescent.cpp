#include "gradientDescent.h"

/*!
move every particle in the direction of the negative energy gradient, scaled by delatT
*/
void gradientDescent::performUpdate()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
    sim->computeForces();
    for (int ii = 0; ii < Ndof; ++ii)
        {
	//cout << "particle " << ii << endl; //debug statement to make sure force calculations are finishing
        displacements[ii] = deltaT*model->forces[ii];
        //printf("p %i d (%f %f %f)  %f\n",ii,displacements[ii][0],displacements[ii][1],displacements[ii][2], displacements[ii].squared_length());
        }
    sim->moveParticles(displacements);
    };
