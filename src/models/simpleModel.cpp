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

    neighborStructure = make_shared<baseNeighborStructure>();
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
    neighborStructure->setInteractionRange(maximumInteractionRange);
    neighborStructure->initialize(positions);
    for (int ii =0; ii < N; ++ii)
        {
        //use the neighborStructure to find candidate neighbors. This function will fill the indices of the neighbors[ii] data structure and populate the positions of the corresponding targetParticles
        vector<meshPosition> targetParticles;
        neighborStructure->constructCandidateNeighborList(positions[ii], ii, neighbors[ii], targetParticles);

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
