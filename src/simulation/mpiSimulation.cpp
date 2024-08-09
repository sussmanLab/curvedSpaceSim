#include "mpiSimulation.h"
/*! \file mpiSimulation.cpp */

/*!
This function sets up MPI_Barriers that ask the mpiModel associated with the
simulation to populate a buffer of integers and a buffer of doubles to be
globally communicated.  For instance, this might be particle type and position
information, but it is up to the model to decide what needs to be shared.  MPI
routines then synchronize this information across all ranks
*/
void mpiSimulation::synchronizeAndTransferBuffers()
    {
    MPI_Barrier(MPI_COMM_WORLD);
    if(totalRanks > 1)
        {
        auto Conf = mConfiguration.lock();
        Conf->processSendingBuffer();
        int n = Conf->largestNumberOfParticlesPerRank;
//printf("r %i, size %i %i\n", myRank,Conf->intTransferBufferSend.size(),Conf->intTransferBufferReceive.size());
        MPI_Allgather(&(Conf->intTransferBufferSend)[0],n,MPI_INT,
                      &(Conf->intTransferBufferReceive)[0],n,MPI_INT,
                      MPI_COMM_WORLD);
        MPI_Allgather(&(Conf->doubleTransferBufferSend[0]),3*n,MPI_DOUBLE,
                      &(Conf->doubleTransferBufferReceive[0]),3*n,MPI_DOUBLE,
                      MPI_COMM_WORLD);
        Conf->processReceivingBuffer();
        }
    else //correctly handle case of mpirun -np 1
        {
        auto Conf = mConfiguration.lock();
        Conf->processSendingBuffer();
        int n = Conf->largestNumberOfParticlesPerRank;
        for(int ii = 0; ii < n; ++ii)
            {
            Conf->intTransferBufferReceive[ii] = Conf->intTransferBufferSend[ii];
            Conf->doubleTransferBufferReceive[3*ii+0] = Conf->doubleTransferBufferSend[3*ii+0];
            Conf->doubleTransferBufferReceive[3*ii+1] = Conf->doubleTransferBufferSend[3*ii+1];
            Conf->doubleTransferBufferReceive[3*ii+2] = Conf->doubleTransferBufferSend[3*ii+2];
            };
        Conf->processReceivingBuffer();
        }
    transfersUpToDate = true;
    }

/*!
Has all models move the particles they are responsible for, then places a call
to synchronize any information needed.
*/
void mpiSimulation::moveParticles(vector<vector3> &displacements)
    {
        {
    auto Conf = mConfiguration.lock();
    Conf->moveParticles(displacements);
        }
    transfersUpToDate = false;
    synchronizeAndTransferBuffers();
    };

/*!
Set a pointer to the configuration
*/
void mpiSimulation::setConfiguration(MPIConfigPtr _config)
    {
    mConfiguration = _config;
    };

/*!
For now: in your main cpp file make a database class on rank 0, and use it (and
its global information) to save the state.
*/
void mpiSimulation::saveState(string fname)
    {
    UNWRITTENCODE("saving states in a mpi simulation not written");
    };

void mpiSimulation::manipulateUpdaterData(vector<double> &data, function<double(double, double)> manipulatingFunction)
    {
    if(totalRanks >1)
        {
        int elements = data.size();
        int rElements = elements * totalRanks;
        dataBuffer = vector<double>(rElements,0.0);
        MPI_Allgather(&data[0],elements,MPI_DOUBLE,&dataBuffer[0],elements,MPI_DOUBLE,MPI_COMM_WORLD);
        
        for (int ii = 0; ii < elements; ++ii) 
            data[ii] = 0.0;

        for (int ii = 0; ii < elements; ++ii)
            {
            for (int rr = 0; rr < totalRanks; ++rr)
                {
                data[ii] = manipulatingFunction(data[ii],dataBuffer[rr*elements+ii]);
                }
            };
        };
    };


