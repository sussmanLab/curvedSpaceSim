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
    neighbors.resize(n);
    neighborVectors.resize(n);
    neighborDistances.resize(n);
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

void simpleModel::findNeighbors(double maximumInteractionRange)
    {
    //pre-optimization, just return an all-to-all coupling matrix
    for (int ii =0; ii < N; ++ii)
        {
        vector<int> currentNeighborList; 
        vector<meshPosition> targetParticles;
        for(int jj = 0; jj < N; ++jj)
            {
            if(ii!=jj)
                {
                currentNeighborList.push_back(jj);
                targetParticles.push_back(positions[jj]);
                };
            }
        //populate the list of neighbor indexes
        neighbors[ii] = currentNeighborList;
        //additionally use the space to populate the list of distances and separation vectors
        vector<vector3> placeholderVector;
        vector<vector3> tangentVector;
        vector<double> distances;
        space->distance(positions[ii],targetParticles,distances,tangentVector,placeholderVector);
        neighborDistances[ii] = distances;
        neighborVectors[ii] = tangentVector;
        };
    };


void simpleModel::setParticlePositions(vector<meshPosition> &newPositions)
    {
    if(N !=newPositions.size())
        initializeSimpleModel(newPositions.size());
    for (int pp = 0;pp < N; ++pp)
        {
        positions[pp] = newPositions[pp];
        };
    };


void simpleModel::computeForces(bool zeroOutForces)
    {
    if(zeroOutForces)
        {
        for(int ii = 0; ii < N; ++ii)
            {
            forces[ii] = vector3(0.0,0.0,0.0);
            };
        }
    };
