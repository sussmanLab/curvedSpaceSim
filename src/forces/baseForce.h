#ifndef baseForce_H
#define baseForce_H

#include "std_include.h"
#include "debuggingHelp.h"
#include "simpleModel.h"
#include "basicSimulation.h"

/*! \file baseForce.h */
//!A base class for implementing force calculations
/*!
In this package we currently assume that we are working with mostly pairwise
interactions, derived classes just need to implement the pairwiseEnergy and
pairwiseForce functions.
*/
class force
    {
    public:
        force(){};
        virtual ~force() = default;

        //!implement the ability to report the name of the class
        virtual string reportSelfName(){string ans = "base force"; return ans;};

        //!the call to compute forces, and store them in the referenced variable
        virtual void computeForces(vector<vector3> &forces,bool zeroOutForce = true, int type = 0);        
        //!some generic function to set parameters...it is the responsibility of child classes to implement this correctly
        virtual void setForceParameters(vector<double> &params){};

        //! A pointer to the governing simulation
        SimPtr sim;
        //!Associate the simulation with this class
        void setSimulation(shared_ptr<basicSimulation> _sim){sim=_sim;};

        //! virtual function to allow the associated model to be a derived class
        virtual void setModel(shared_ptr<simpleModel> _model){model=_model;};

        //!compute the energy associated with this force
        virtual double computeEnergy(bool verbose = false);
        //! child classes must implement a specific function that, given a unit vector specifying the direction of separation and scalar distance, returns an energy
        virtual double pairwiseEnergy(vector3 separation,double distance) = 0;
        //! child classes must implement a specific function that, given a unit vector specifying the direction of separation and scalar distance, returns a vector force
        virtual vector3 pairwiseForce(vector3 separation,double distance) = 0;

        //! A pointer to a simpleModel that the updater acts on
        shared_ptr<simpleModel> model;

        //!Forces might update the total energy associated with them
        double energy;
        //!forces might have a maximum range of interaction
        double maximumInteractionRange=1;

        //!TODO: update getClassSizes for all class
        virtual double getClassSize()
            {
            return 0.000000001*((2+0)*sizeof(double));
            }
    };

typedef shared_ptr<force> ForcePtr;
typedef weak_ptr<force> WeakForcePtr;
#endif
