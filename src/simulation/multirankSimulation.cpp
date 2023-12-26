#include "mpiSimulation.h"
/*! \file mpiSimulation.cpp */

void mpiSimulation::sumUpdaterData(vector<scalar> &data)
    {
    /*
    if(nRanks >1)
        {
        int elements = data.size();
        int rElements = elements * nRanks;
        if (dataBuffer.size() < rElements)
            dataBuffer.resize(rElements);
        p1.start();
        MPI_Allgather(&data[0],elements,MPI_SCALAR,&dataBuffer[0],elements,MPI_SCALAR,MPI_COMM_WORLD);
        p1.end();
        for (int ii = 0; ii < elements; ++ii) data[ii] = 0.0;

        for (int ii = 0; ii < elements; ++ii)
            for (int rr = 0; rr < nRanks; ++rr)
                {
                data[ii] += dataBuffer[rr*elements+ii];
                }
        };
    */
    };

void mpiSimulation::synchronizeAndTransferBuffers()
    {
    if(nRanks > 1)
        {
        auto Conf = mConfiguration.lock();
        Conf->processSendingBuffer();
        int n = Conf->N;
        MPI_Allgather(Conf->&intTransferBufferSend[0],n,MPI_INT,
                      Conf->&intTransferBufferReceive[0],n,MPI_INT,
                      MPI_COMM_WORLD);
        MPI_Allgather(Conf->&doubleTransferBufferSend[0],3*n,MPI_DOUBLE,
                      Conf->&doubleTransferBufferReceive[0],3*n,MPI_DOUBLE,
                      MPI_COMM_WORLD);
        Conf->processReceivingBuffer();
        };
        transfersUpToDate = true;
        };
    }

/*!
Calls the configuration to displace the degrees of freedom, and communicates halo sites according
to the rankTopology and boolean settings
*/
void mpiSimulation::moveParticles(GPUArray<dVec> &displacements,scalar scale)
    {
        {
    auto Conf = mConfiguration.lock();
    Conf->moveParticles(displacements,scale);
        }
    transfersUpToDate = false;
    synchronizeAndTransferBuffers();
    };

/*!
Set a pointer to the configuration
*/
void mpiSimulation::setConfiguration(MConfigPtr _config)
    {
    mConfiguration = _config;
    };

/*!
*/
void mpiSimulation::saveState(string fname)
    {
    UNWRITTENCODE("saving states in a mpi simulation not written");
    };
