#ifndef MPIMODEL_H
#define MPIMODEL_H

#include "simpleModel.h"
#include <mpi.h>

/*! \file mpiModel.h
 */

/*!
The simplest implementation of an MPI extension of the simpleModel class... rank 0 is used to set any global information (say, initial particle positions), and then each rank gets approximately (N/(numberOfRanks)) particles to take care of. Coordinates tightly with the functions in mpiSimulation.

In this implementation, only positional information is shared.  A separate class
would be needed if velocity information also is to be coordinated across ranks
(as in, e.g., an updater in which velocities locally align).
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

        //!On mpiModels, this uses rank 0 to *globally* set all particle positions
        virtual void setRandomParticlePositions(noiseSource &noise);
        //!On mpiModels, this uses rank 0 to *globally* set all particle velocities
        virtual void setMaxwellBoltzmannVelocities(noiseSource &noise, double T);

        //!have (by default) rank 0 fill and then broadcast the globalPositions vector
        void broadcastParticlePositions(vector<meshPosition> &p, int broadcastRoot = 0);
        //!have (by default) rank 0 fill and then broadcast the globalVelocities vector... mostly useful for setting initial conditions
        void broadcastParticleVelocities(vector<vector3> &v, int broadcastRoot = 0);
        
        //!before MPI communication, prepare the sending buffer with relevant model data
        virtual void processSendingBuffer(int directionType = -1); 
        //! after MPI communication, convert the receiving buffer into the relevant model data
        virtual void processReceivingBuffer(int directionType = -1);  

        //!look at the global position information and update the information that should be local to this process
        void readFromGlobalPositions();

        //!specialize the findNeighbors routine to deal with particle index offsetting
        virtual void findNeighbors(double maximumInteractionRange);

        //!The total number of particles across all ranks
        int NTotal;

        //!global particle  positions (i.e., the positions of *all* particles in the simulation)
        vector<meshPosition> globalPositions;
        //!global particle velocities
        vector<vector3> globalVelocities;

        //Buffers for sending and receiving ints and doubles for inter-process communication
        vector<int> intTransferBufferSend;
        vector<double> doubleTransferBufferSend;
        vector<int> intTransferBufferReceive;
        vector<double> doubleTransferBufferReceive;

        int largestNumberOfParticlesPerRank;
        //!convert from the internal meshPosition representation and fill the vector of euclidean locations
        virtual void fillEuclideanLocations();

    protected:
        //the current rank has N of NTotal degrees of freedom in the globalPositions vector, corresponding to "for (int ii = minIdx; ii < maxIdx; ++ii)" type accesses
        //!the global index corresponding to the current rank's "particle 0"
        int minIdx;
        //!the global index corresponding to the current rank's last particle
        int maxIdx;
        //!The local rank of the process
        int localRank;
        //! the total number of processes that are communicating in the simulation
        int totalRanks;

        //!Use nTotal, localRank, and totalRanks to set min/maxIdx and N
        virtual void determineIndexBounds();

    };
typedef shared_ptr<mpiModel> MPIConfigPtr;
typedef weak_ptr<mpiModel> MPIWeakConfigPtr;
#endif
