#ifndef updater_H
#define updater_H

#include "std_include.h"
#include "simpleModel.h"
#include "basicSimulation.h"

/*! \file baseUpdater.h */
//!A base class for implementing simple updaters
/*!
An updater is some class object that can update something about the underlying state of the system. An example might be an equation of motion, or an updater that periodically subtracts off
any center-of-mass motion of a system as it evolves, etc.. A simulation will call all updaters in a loop, e.g. for(each updater i in list) updater[i].Update(Timestep) To facilitate this structure, but acknowledge that any given updater might only need to be called occasionally, the Update function is passed a timestep, and each updaters has a period that should be set.
*/
class updater
    {
    public:
        //! by default, updaters are called every timestep with no offset
        updater(){Period = -1;Phase = 0;reproducible = true;};
        //! construct an updater with a specific period
        updater(int _p){Period = _p; Phase = 0;};
        virtual ~updater() = default;
        //! The fundamental function that a controlling Simulation can call
        virtual void Update(int timestep);
        //! The function which performs the update
        virtual void performUpdate();
        //! A pointer to the governing simulation
        shared_ptr<basicSimulation> sim;
        //!set the simulation
        void setSimulation(shared_ptr<basicSimulation> _sim){sim=_sim;};

        //! A pointer to a simpleModel that the updater acts on
        shared_ptr<simpleModel> model;
        //! virtual function to allow the model to be a derived class
        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            initializeFromModel();
            };

        //!by default, set Ndof from the number of particles in the model
        virtual void initializeFromModel(){Ndof = model->getNumberOfParticles();};

        //! set the period
        void setPeriod(int _p){Period = _p;};
        //! set the phase
        void setPhase(int _p){Phase = _p;};

        //!allow for spatial sorting to be called if necessary...
        virtual void spatialSorting(){};

        //!Allow for a reproducibility call to be made
        virtual void setReproducible(bool rep){reproducible = rep;};

        //!Get the number of degrees of freedom of the equation of motion
        int getNdof(){return Ndof;};
        //!Set the number of degrees of freedom of the equation of motion
        void setNdof(int _n){Ndof = _n;};

        //!allow all updaters to potentially implement an internal time scale
        virtual void setDeltaT(double dt){deltaT = dt;};

        //!allow for setting multiple threads
        virtual void setNThreads(int n){nThreads = n;};
        //!A sample function, mostly to show how to use manipulateUpdaterData
        virtual double getMaxForce();
        //!A sample function, mostly to show how to use manipulateUpdaterData
        virtual double getForceNorm();

        //!Set the maximum number of iterations before terminating
        void setMaximumIterations(int maxIt=-1){maxIterations = maxIt;};
        int getCurrentIterations(){return iterations;};
        int getMaxIterations(){return maxIterations;};
        void setCurrentIterations(int newIterations){iterations=newIterations;};
        //! TODO
        virtual double getClassSize()
            {
            return 0.000000001*(6*sizeof(int) + 2*sizeof(bool) + sizeof(double));
            }

        //!The number of iterations performed
        int iterations;

        double squaredTotalForceNorm;
        double maximumForceNorm;

    protected:
        //!number of threads to use TODO
        int nThreads=1;
        //!The period of the updater... the updater will work every Period timesteps
        int Period;
        //!The phase of the updater... the updater will work every Period timesteps offset by a phase
        int Phase;
        //!some measure of the number of degrees of freedom the equations of motion might need to know about locally
        int Ndof;
        //!whether the RNGs give reproducible results
        bool reproducible;
        //!The internal time step size
        double deltaT;
        //!The maximum number of iterations allowed...perhaps useful for energy minimizers, for instance
        int maxIterations = -1;
    };

typedef shared_ptr<updater> UpdaterPtr;
typedef weak_ptr<updater> WeakUpdaterPtr;
#endif
