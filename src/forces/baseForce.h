#ifndef baseForce_H
#define baseForce_H

#include "std_include.h"
#include "simpleModel.h"
#include "basicSimulation.h"

/*! \file baseForce.h */
//!A base class for implementing force calculations
/*!
 *
 *
*/
class force
    {
    public:
        force();

        virtual string reportSelfName(){string ans = "unnamed"; return ans;};

        //!the call to compute forces, and store them in the referenced variable
        virtual void computeForces(vector<vector3> &forces,bool zeroOutForce = true, int type = 0);        
        //!some generic function to set parameters
        virtual void setForceParameters(vector<double> &params);

        //! A pointer to the governing simulation
        SimPtr sim;
        //!set the simulation
        void setSimulation(shared_ptr<basicSimulation> _sim){sim=_sim;};

        //! virtual function to allow the model to be a derived class
        virtual void setModel(shared_ptr<simpleModel> _model){model=_model;};

        //!compute the energy associated with this force
        virtual double computeEnergy(bool verbose = false){return 0.;};

        //! A pointer to a simpleModel that the updater acts on
        shared_ptr<simpleModel> model;

        //!Forces might update the total energy associated with them
        double energy;

        virtual double getClassSize()
            {
            return 0.000000001*(2*sizeof(bool) + (1+0)*sizeof(double));
            }
    };

typedef shared_ptr<force> ForcePtr;
typedef weak_ptr<force> WeakForcePtr;
#endif
