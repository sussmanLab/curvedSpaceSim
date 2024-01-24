#ifndef MPIMODEL_H
#define MPIMODEL_H

#include "simpleModel.h"
#include <mpi.h>

/*! \file mpiModel.h
 */

/*!
*/
class mpiModel : public simpleModel
    {
    public:
        //!The base constructor requires the total number of particles and information about how many ranks there are
        mpiModel(int nTotal, int _localRank, int _totalRanks, bool verbose = false);
        //!a blank default constructor
        mpiModel(){};
        //!initialize the size of the basic data structure arrays
        void initializeMPIModel(int n, int nTotal);

        //!if an mpiModel is created with the blank constructor, be sure to setRankInformation before calling initializeMPIModel
        void setRankInformation(int nTotal, int _localRank, int _totalRanks)
            {
            NTotal = nTotal;
            localRank=_localRank;
            totalRanks=_totalRanks;
            determineIndexBounds();
            };

        //!have (by default) rank 0 fill and then broadcast the globalPositions vector
        void broadcastParticlePositions(vector<meshPosition> &p, int broadcastRoot = 0);
        //!have (by default) rank 0 fill and then broadcast the globalVelocities vector... mostly useful for setting initial conditions
        void broadcastParticleVelocities(vector<vector3> &v, int broadcastRoot = 0);

        virtual void processSendingBuffer(int directionType = -1); 
        virtual void processReceivingBuffer(int directionType = -1);  

        void readFromGlobalPositions();

        virtual void findNeighbors(double maximumInteractionRange);

        //!The total number of particles across all ranks
        int NTotal;

        //!global particle  positions
        vector<meshPosition> globalPositions;

        vector<int> intTransferBufferSend;
        vector<double> doubleTransferBufferSend;
        vector<int> intTransferBufferReceive;
        vector<double> doubleTransferBufferReceive;

        int largestNumberOfParticlesPerRank;
        virtual void fillEuclideanLocations();

    protected:
        //the current rank has N of NTotal degrees of freedom in the globalPositions vector, corresponding to "for (int ii = minIdx; ii < maxIdx; ++ii)" type accesses
        int minIdx;
        int maxIdx;
        int localRank;
        int totalRanks;

        //!Use nTotal, localRank, and totalRanks to set min/maxIdx
        virtual void determineIndexBounds();

        vector<vector3> globalVelocities;
    };
typedef shared_ptr<mpiModel> MPIConfigPtr;
typedef weak_ptr<mpiModel> MPIWeakConfigPtr;
#endif
