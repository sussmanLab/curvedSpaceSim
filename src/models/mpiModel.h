#ifndef MPIMODEL_H
#define MPIMODEL_H

#include "simpleModel.h"

/*! \file mpiModel.h
 */

/*!
*/
class mpiModel : public simpleModel
    {
    public:
        //!The base constructor requires the total number of particles and information about how many ranks there are
        mpiModel(int nTotal, int _localRank, int _totalRanks);
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

        virtual void processSendingBuffer(int directionType = -1); 
        virtual void processReceivingBuffer(int directionType = -1);  


        virtual void findNeighbors(double maximumInteractionRange);

        //!The total number of particles across all ranks
        int NTotal;

        //!global particle  positions
        vector<meshPosition> globalPositions;

        vector<int> intTransferBufferSend;
        vector<double> doubleTransferBufferSend;
        vector<int> intTransferBufferReceive;
        vector<double> doubleTransferBufferReceive;

    protected:
        //the current rank has N of NTotal degrees of freedom in the globalPositions vector, corresponding to "for (int ii = minIdx; ii < maxIdx; ++ii)" type accesses
        int minIdx;
        int maxIdx;
        int localRank;
        int totalRanks;
        int largestNumberOfParticlesPerRank;

        //!Use nTotal, localRank, and totalRanks to set min/maxIdx
        virtual void determineIndexBounds();

    };
typedef shared_ptr<mpiModel> MPIConfigPtr;
typedef weak_ptr<mpiModel> MPIWeakConfigPtr;
#endif
