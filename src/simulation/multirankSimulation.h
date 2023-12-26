#ifndef multirankSIMULATION_H
#define multirankSIMULATION_H

#include "simulation.h"
#include "mpiModel.h"
#include <mpi.h>

/*! \file multirankSimulation.h */
/*!
At the moment: this simple all-to-all multirank simulation framework assumes that all 
"mpiModels" want an all-to-all communcation pattern, and appropriately handle
processSendingBuffer() and processReceivingBuffer(), filling and/or extracting information
from intTransferBufferSend, doubleTransferBufferSend, intTransferBufferReceive, doubleTransferBufferReceive

*/
class multirankSimulation : public simulation, public enable_shared_from_this<multirankSimulation>
    {
    public:
        multirankSimulation(int _myRank,int _totalRanks)
            {
            myRank = _myRank;
            totalRanks=_totalRanks;
            }

        //! synchronize mpi and make transfer buffers
        virtual void synchronizeAndTransferBuffers();

        //!The configuration of latticeSites
        weak_ptr<mpiModel> mConfiguration;

        //!return a shared pointer to this Simulation
        shared_ptr<multirankSimulation> getPointer(){ return shared_from_this();};
        //!Pass in a reference to the configuration
        void setConfiguration(shared_ptr<mpiModel> _config);

        //! manipulate data from updaters
        virtual void sumUpdaterData(vector<scalar> &data);

        void saveState(string fname);

        virtual void reportSelf(){cout << "in the multirank simulation class" << endl;};

    protected:
        int myRank;
        int totalRanks;
        MPI_Status mpiStatus;
        vector<MPI_Status> mpiStatuses;
        vector<MPI_Request> mpiRequests;
        vector<scalar> dataBuffer;
        bool transfersUpToDate;
    };
#endif
