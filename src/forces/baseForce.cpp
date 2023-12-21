#include "baseForce.h"
/*! \file baseForce.cpp */

force::force()
    {
    };

void force::setForceParameters(vector<double> &params)
    {
    };

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
            forces[ii] += pairwiseForce(model->neighborVectors[ii][jj],model->neighborDistances[ii][jj]);
        }
    };

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
    };
