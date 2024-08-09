#ifndef mpiSIMULATION_H
#define mpiSIMULATION_H

#include "simulation.h"
#include "mpiModel.h"
#include <mpi.h>

/*! \file mpiSimulation.h */

//! A simulation class that evenly splits up all particles between mpi ranks (no spatial sorting)
/*!
At the moment: this simple all-to-all mpi simulation framework assumes that all
"mpiModels" want an all-to-all communcation pattern.  The simple strategy here
is that all ranks maintain a global list of all particle positions; each rank is
then responsible for computing forces only for a subset of these particles.  Any
time particles are moved, mpi communication happens (via processSendingBuffer()
and processReceivingBuffer(), filling and/or extracting information from
intTransferBufferSend, doubleTransferBufferSend, intTransferBufferReceive,
doubleTransferBufferReceive) to make sure the set of particle information known
to the mpi rank is up-to-date.  This pattern is obviously inefficient for
normal, flat-space simulations of finite-range forces.  However, in
curvedSpaceSim so much of the computation time is spent computing geodesic paths
that splitting up just this task among ranks is a reasonable first choice.
*/
class mpiSimulation : public simulation
    {
    public:
        mpiSimulation(int _myRank,int _totalRanks)
            {
            myRank = _myRank;
            totalRanks=_totalRanks;
            }
        //!Most importantly, moveParticles is modified to also invoke MPI communication between ranks
        virtual void moveParticles(vector<vector3> &displacements);

        //! synchronize mpi and make transfer buffers
        virtual void synchronizeAndTransferBuffers();

        shared_ptr<mpiSimulation> getPointer()
            {
            return dynamic_pointer_cast<mpiSimulation>(simulation::shared_from_this());
            }

        //!The configuration
        weak_ptr<mpiModel> mConfiguration;

        //!Pass in a reference to the configuration
        void setConfiguration(shared_ptr<mpiModel> _config);

        //! manipulate data from updaters, communicating this information to all ranks
        virtual void manipulateUpdaterData(vector<double> &data, function<double(double, double)> manipulatingFunction);

        //!TODO: saving MPI state data is currently not implemented through the simulation class...
        void saveState(string fname);

        virtual void reportSelf(){cout << "in the mpi simulation class" << endl;};

    protected:
        MPI_Status mpiStatus;
        vector<MPI_Status> mpiStatuses;
        vector<MPI_Request> mpiRequests;
        vector<double> dataBuffer;
        bool transfersUpToDate;
    };
#endif
