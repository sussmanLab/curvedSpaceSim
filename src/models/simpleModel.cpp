#include "simpleModel.h"
/*! \file simpleModel.cpp */

/*!
 * Set the size of basic data structures...
*/
simpleModel::simpleModel(int n) :
    N(n)
    {
    cout << "initializing a model with "<< N << " particles" << endl;
    initializeSimpleModel(n);
    };

/*!
 * actually set the array sizes. positions, velocities, forces are zero
 * masses are set to unity
*/
void simpleModel::initializeSimpleModel(int n)
    {
    N=n;

    positions.resize(n);
    velocities.resize(n);
    forces.resize(n);
    //types.resize(n);
    //masses.resize(n);
    //radii.resize(n);
    //vector<scalar> halves(N,.5);
    //vector<int> units(N,0);
    };

void simpleModel::moveParticles(vector<vector3> &disp)
    {
    for(int ii = 0; ii < N; ++ii)
        {
        space->displaceParticle(positions[ii],disp[ii]);
        };
    };

/*
void simpleModel::setParticlePositions(vector<dVec> &newPositions)
    {
    if(N !=newPositions.size())
        initializeSimpleModel(newPositions.size());
    ArrayHandle<dVec> p(positions);
    for (int pp = 0;pp < N; ++pp)
        {
        p.data[pp] = newPositions[pp];
        Box->putInBoxReal(p.data[pp]);
        };
    };
*/

void simpleModel::computeForces(bool zeroOutForces)
    {
    };
