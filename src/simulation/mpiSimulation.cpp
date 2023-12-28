#include "mpiSimulation.h"
/*! \file mpiSimulation.cpp */

void mpiSimulation::sumUpdaterData(vector<double> &data)
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
        MPI_Allgather(&(Conf->intTransferBufferSend)[0],n,MPI_INT,
                      &(Conf->intTransferBufferReceive)[0],n,MPI_INT,
                      MPI_COMM_WORLD);
        MPI_Allgather(&(Conf->doubleTransferBufferSend[0]),3*n,MPI_DOUBLE,
                      &(Conf->doubleTransferBufferReceive[0]),3*n,MPI_DOUBLE,
                      MPI_COMM_WORLD);
        Conf->processReceivingBuffer();
        };
    transfersUpToDate = true;
    }

void mpiSimulation::computeForces()
    {
    auto Conf = mConfiguration.lock();
    for (unsigned int f = 0; f < forceComputers.size(); ++f)
        {
        auto frc = forceComputers[f].lock();
        bool zeroForces = (f==0);
        frc->computeForces(Conf->forces,zeroForces);
        };
    };
/*!
Calls the configuration to displace the degrees of freedom, and communicates halo sites according
to the rankTopology and boolean settings
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
*/
void mpiSimulation::saveState(string fname)
    {
    UNWRITTENCODE("saving states in a mpi simulation not written");
    };
