#include "mpiModel.h"
/*! \file mpiModel.cpp */

/*!
 * Set the size of basic data structures...
*/
mpiModel::mpiModel(int nTotal, int _localRank, int _totalRanks) :
    NTotal(nTotal), localRank(_localRank), totalRanks(_totalRanks) 
    {
    determineIndexBounds();
    N = maxIdx-minIdx;
    cout << "initializing a model with "<< N << " particles" << endl;
    initializeMPIModel(N,NTotal);
    };

void mpiModel::determineIndexBounds()
    {
    largestNumberOfParticlesPerRank = ceil((double) nTotal / (double) totalRanks);
    minIdx = localRank*largestNumberOfParticlesPerRank;
    maxIdx = (localRank+1)*largestNumberOfParticlesPerRank;
    //The last rank handles a few fewer particles due to rounding
    if(localRank = totalRanks-1)
        maxIdx = nTotal;
    N = maxIdx-minIdx;
    }
/*!
 * actually set the array sizes. positions, velocities, forces are zero
 * masses are set to unity
*/
void mpiModel::initializeMPIModel(int n, int nTotal)
    {
    initializeSimpleModel(n);

    NTotal = nTotal;
    globalPositions.resize(NTotal);

    //transfer buffers are padded to all be the same size
    intTransferBufferSend.resize(largestNumberOfParticlesPerRank);
    doubleTransferBufferSend.resize(3*largestNumberOfParticlesPerRank);

    intTransferBufferReceive.resize(largestNumberOfParticlesPerRank*totalRanks);
    doubleTransferBufferReceive.resize(3*largestNumberOfParticlesPerRank*totalRanks);
    };

void mpiModel::findNeighbors(double maximumInteractionRange)
    {
    //pre-optimization, just return an all-to-all coupling matrix

    //i are particles on this rank; j are particles in the global vector
    for (int ii =0; ii < N; ++ii)
        {
        vector<int> currentNeighborList; 
        vector<meshPosition> targetParticles;
        for(int jj = 0; jj < NTotal; ++jj)
            {
            if(ii!=jj)
                {
                currentNeighborList.push_back(jj);
                targetParticles.push_back(globalPositions[jj]);
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

void mpiModel::processSendingBuffer(int directionType)
    {
    for (int ii = 0; ii < N; ++ii)
        {
        intTransferBufferSend[ii] = positions[ii].faceIndex;
        doubleTransferBufferSend[3*ii+0] = positions[ii].x[0];
        doubleTransferBufferSend[3*ii+1] = positions[ii].x[1];
        doubleTransferBufferSend[3*ii+2] = positions[ii].x[2];
        };
    };

void mpiModel::processReceivingBuffer(int directionType)
    {
    double x,y,z;
    for (int ii = 0; ii < NTotal; ++ii)
        {
        x = doubleTransferBufferReceive[3*ii+0];
        y = doubleTransferBufferReceive[3*ii+1];
        z = doubleTransferBufferReceive[3*ii+2];
        globalPositions[ii].x = point3(x,y,z);
        globalPoitions.faceIndex = intTransferBufferReceive[ii];
        };
    for (int ii = 0; ii < N; ++ii)
        {
        positions[ii].x = point3(globalPositions[ii+minIdx][0],globalPositions[ii+minIdx][1],globalPositions[ii+minIdx][2]);
        positions[ii].faceIndex=globalPositions[ii].faceIndex;
        }
    };

