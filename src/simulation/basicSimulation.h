#ifndef BASICSIMULATION_H
#define BASICSIMULATION_H
/*! \file basicSimulation.h */
//!Basic simulations just know that there are virtual functions that implement computeForces, moveParticles, and have a simulation domain

#include "simpleModel.h"
#include <functional>

class basicSimulation
    {
    public:
        //!Initialize all the shared pointers, etc.
        basicSimulation();

        //!Call the force computer to compute the forces
        virtual void computeForces()=0;
        //!Call the configuration to move particles around
        virtual void moveParticles(vector<vector3> &displacements)=0;
        //!compute the potential energy associated with all of the forces
        virtual double computePotentialEnergy(bool verbose =false){return 0.0;};

        //!The configuration of particles
        WeakConfigPtr configuration;
        //! An integer that keeps track of how often performTimestep has been called
        int integerTimestep;
        //!The current simulation time
        double Time;
        //! The dt of a time step
        double integrationTimestep;

        //!Set the time between spatial sorting operations.
        void setSortPeriod(int sp){sortPeriod = sp;};

        //!reset the simulation clock
        virtual void setCurrentTime(double _cTime);
        //!reset the simulation clock counter
        virtual void setCurrentTimestep(int _cTime){integerTimestep =_cTime;};

        //! manipulate data from updaters...on multirank simulations, might be needed to get, e.g., current total KE, or the current size of the force vector
        virtual void manipulateUpdaterData(vector<double> &data, function<double(double, double)> manipulatingFunction){};
        //!integer for this rank (if simulations are being done in parallel
        int myRank;
        //!total number of ranks
        int totalRanks;

        virtual void reportSelf(){cout << "in the base simulation class" << endl;cout.flush();};
    protected:
        int sortPeriod;
    };
typedef shared_ptr<basicSimulation> SimPtr;
typedef weak_ptr<basicSimulation> WeakSimPtr;

#endif
