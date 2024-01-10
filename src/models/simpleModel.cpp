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
    euclideanLocations.resize(n);
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

void simpleModel::fillEuclideanLocations()
    {
    space->meshPositionToEuclideanLocation(positions,euclideanLocations);
    }

void simpleModel::moveParticles(vector<vector3> &disp)
    {
    for(int ii = 0; ii < N; ++ii)
        {
        space->displaceParticle(positions[ii],disp[ii]);
        };
    };

void simpleModel::findNeighbors(double maximumInteractionRange)
    {
    //first, determine any needed neighbor structure initialization
    neighborStructure->setInteractionRange(maximumInteractionRange);
    bool euclideanNeighborsMeshPositions = neighborStructure->requireEuclideanPositions && !space->positionsAreEuclidean;
    if(euclideanNeighborsMeshPositions)
        {
        space->meshPositionToEuclideanLocation(positions,euclideanMeshPosition);
        neighborStructure->initialize(euclideanMeshPosition);
        }
    else
        neighborStructure->initialize(positions);

    double meanNN = 0;
    for (int ii =0; ii < N; ++ii)
        {
        //use the neighborStructure to find candidate neighbors. This function will fill the indices of the neighbors[ii] data structure and populate the positions of the corresponding targetParticles
        vector<meshPosition> targetParticles;
        //if needed, target a euclidean position
        if(euclideanNeighborsMeshPositions)
            neighborStructure->constructCandidateNeighborList(euclideanMeshPosition[ii], ii, neighbors[ii], targetParticles);
        else
            neighborStructure->constructCandidateNeighborList(positions[ii], ii, neighbors[ii], targetParticles);

        //if needed, re-fill targetParticles with space-appropriate data
        if(euclideanNeighborsMeshPositions)
            {
            for(int jj = 0; jj < neighbors[ii].size();++jj)
                targetParticles[jj] = positions[neighbors[ii][jj]];
            };
        //additionally use the space to populate the list of distances and separation vectors
        vector<vector3> placeholderVector;
        vector<vector3> tangentVector;
        vector<double> distances;
        if(neighborStructure->requireEuclideanPositions && !(space->positionsAreEuclidean))
            {
            space->distance(positions[ii],targetParticles,distances,tangentVector,placeholderVector);
            }
        else
            space->distance(positions[ii],targetParticles,distances,tangentVector,placeholderVector);
        neighborDistances[ii] = distances;
        neighborVectors[ii] = tangentVector;
        meanNN += distances.size();
        };
    if(verbose)
        printf("mean number of neighbors = %f\n",meanNN/N);
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
