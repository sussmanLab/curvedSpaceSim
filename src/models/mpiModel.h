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
        //!The base constructor requires the number of particles
        mpiModel(int n, int nTotal);
        //!a blank default constructor
        mpiModel(){};
        //!initialize the size of the basic data structure arrays
        void initializeMPIModel(int n, int nTotal);

        virtual void prepareSendingBuffer(int directionType = -1); 
        virtual void readReceivingBuffer(int directionType = -1);  


        virtual void findNeighbors(double maximumInteractionRange);

        //!The total number of particles across all ranks
        int NTotal;

        //!global particle  positions
        vector<meshPosition> globalPositions;


    protected:
        shared_ptr<baseSpace> space;
        //the current rank has N of NTotal degrees of freedom in the globalPositions vector, corresponding to "for (int ii = minIdx; ii < maxIdx; ++ii)" type accesses
        int minIdx;
        int maxIdx;
        vector<int> intTransferBufferSend;
        vector<double> doubleTransferBufferSend;
        vector<int> intTransferBufferReceive;
        vector<double> doubleTransferBufferReceive;

    };
typedef shared_ptr<mpiModel> MPIConfigPtr;
typedef weak_ptr<mpiModel> MPIWeakConfigPtr;
#endif
