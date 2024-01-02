#include "mpiModel.h"
/*! \file mpiModel.cpp */

/*!
 * Set the size of basic data structures...
*/
mpiModel::mpiModel(int nTotal, int _localRank, int _totalRanks, bool verbose)
    {
    NTotal = nTotal;
    localRank=_localRank;
    totalRanks = _totalRanks;
    determineIndexBounds();
    if(verbose)
        cout << "initializing a model on rank " << _localRank << " of " << _totalRanks <<"  with "<< N << " particles" << endl;
    initializeMPIModel(N,NTotal);
//printf("r %i, minIdx maxIdx=(%i,%i)\n ",localRank,minIdx,maxIdx);
    };

void mpiModel::determineIndexBounds()
    {
    largestNumberOfParticlesPerRank = ceil((double) NTotal / (double) totalRanks);
    minIdx = localRank*largestNumberOfParticlesPerRank;
    maxIdx = (localRank+1)*largestNumberOfParticlesPerRank;
    //The last rank handles a few fewer particles due to rounding
    if(localRank == totalRanks-1)
        maxIdx = NTotal;
    N = maxIdx-minIdx;
    //printf("rank %i, largestP %i, min %i, max %i, N %i\n",localRank, largestNumberOfParticlesPerRank, minIdx, maxIdx, N);
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
    neighborStructure->setInteractionRange(maximumInteractionRange);
    neighborStructure->initialize(globalPositions);

    //i are particles on this rank; j are particles in the global vector

    for (int ii =0; ii < N; ++ii)
        {
        //use the neighborStructure to find candidate neighbors. This function will fill the indices of the neighbors[ii] data structure and populate the positions of the corresponding targetParticles. minIdx provides the necessary offset for indexing between positions and globalPositions
        vector<meshPosition> targetParticles;
        neighborStructure->constructCandidateNeighborList(positions[ii], ii, neighbors[ii], targetParticles,minIdx);

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
        globalPositions[ii].faceIndex = intTransferBufferReceive[ii];
        };
    readFromGlobalPositions();
    };

void mpiModel::readFromGlobalPositions()
    {
    for (int ii = 0; ii < N; ++ii)
        {
        positions[ii].x = point3(globalPositions[ii+minIdx].x[0],globalPositions[ii+minIdx].x[1],globalPositions[ii+minIdx].x[2]);
        positions[ii].faceIndex=globalPositions[ii+minIdx].faceIndex;
//        printf("r %i p %i (%f,%f,%f)\n",localRank,ii,positions[ii].x[0],positions[ii].x[1],positions[ii].x[2]);
        }
    };

void mpiModel::broadcastParticlePositions(vector<meshPosition> &p, int broadcastRoot)
    {
    if(localRank == broadcastRoot)
        {
        for (int ii = 0; ii < NTotal; ++ii)
            {
            doubleTransferBufferReceive[3*ii+0]=p[ii].x[0];
            doubleTransferBufferReceive[3*ii+1]=p[ii].x[1];
            doubleTransferBufferReceive[3*ii+2]=p[ii].x[2];
            intTransferBufferReceive[ii] = p[ii].faceIndex;
            };
        };
    MPI_Bcast(&intTransferBufferReceive[0],NTotal,MPI_INT,broadcastRoot,MPI_COMM_WORLD);
    MPI_Bcast(&doubleTransferBufferReceive[0],3*NTotal,MPI_DOUBLE,broadcastRoot,MPI_COMM_WORLD);
    processReceivingBuffer();

    };
