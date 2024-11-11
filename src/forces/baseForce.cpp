#include "baseForce.h"
/*! \file baseForce.cpp */

/*!
Assuming that we are implementing pairwise forces, this function will first call
the model function to find neighbors within some maximumInteractionRange.  This
is a vector of vectors containing the index,j, of all neighbors of particle i;
this function executes the double loop and calls the child class to compute the
specific pairwise force between the particles.  Note that computeForces expects
the neighborVectors to all be *unit vectors*
*/
void force::computeForces(vector<vector3> &forces,bool zeroOutForce, int type)
    {
    if(forces.size() != model->N)
        forces.resize(model->N);
    model->findNeighbors(maximumInteractionRange);

    for (int ii = 0; ii < model->N; ++ii)
        {
        if(zeroOutForce)
            forces[ii] = vector3(0.0,0.0,0.0);
        int neighborNumber = model->neighbors[ii].size();
        for (int jj = 0; jj < neighborNumber; ++jj)
            {
            forces[ii] += pairwiseForce(model->neighborVectors[ii][jj],model->neighborDistances[ii][jj]);
            }
        }
    };

/*!
See the comment on computeForces... computeEnergy works the same way
*/
double force::computeEnergy(bool verbose)
    {
    model->findNeighbors(maximumInteractionRange);
    energy = 0;
    for (int ii = 0; ii < model->N; ++ii)
        {
        int neighborNumber = model->neighbors[ii].size();
        for (int jj = 0; jj < neighborNumber; ++jj)
            energy += pairwiseEnergy(model->neighborVectors[ii][jj],model->neighborDistances[ii][jj]);
        }
    return energy;
    };
