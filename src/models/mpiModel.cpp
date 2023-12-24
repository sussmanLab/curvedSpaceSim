#include "mpiModel.h"
/*! \file mpiModel.cpp */

/*!
 * Set the size of basic data structures...
*/
mpiModel::mpiModel(int n, int nTotal) :
    N(n)
    {
    cout << "initializing a model with "<< N << " particles" << endl;
    initializeSimpleModel(n,nTotal);
    };

/*!
 * actually set the array sizes. positions, velocities, forces are zero
 * masses are set to unity
*/
void mpiModel::initializeMPIModel(int n, int nTotal)
    {
    initializeSimpleModel(n);

    NTotal = nTotal;
    globalPositions.resize(NTotal);
    };

void mpiModel::findNeighbors(double maximumInteractionRange)
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
