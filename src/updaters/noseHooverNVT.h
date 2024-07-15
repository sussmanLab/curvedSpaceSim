#ifndef noseHooverNVT_H
#define noseHooverNVT_H

#include "baseUpdater.h"

/*! \file noseHooverNVT.h */
//! Implements NVT dynamics according to the Nose-Hoover equations of motion with a chain of thermostats
/*!
 *This allows one to do standard NVT simulations. A chain (whose length can be specified by the user)
 of thermostats is used to maintain the target temperature. We closely follow the Frenkel & Smit
 update scheme, which is itself based on:
 Martyna, Tuckerman, Tobias, and Klein
 Mol. Phys. 87, 1117 (1996)
*/

class noseHooverNVT : public updater
    {
    public:
        //The constructor allows one to set the time constant associated with the bath variable, but this needs to be implemented in the code by letting tau alter the bath variable mass
        noseHooverNVT(double _dt, double _T, double _tau=1.0, int _M=2);

        virtual void performUpdate();

        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            model->particleShiftsRequireVelocityTransport = true;
            initializeFromModel();
            setBathVariables();
            };
    protected:
        vector<vector3> displacements;
        double temperature;
        int chainLength;
        //!bathVariables contains (pos, vel, accel, mass) of the chain variables
        vector<double4> bathVariables;
        double timeConstant;
        void propagateChain();
        void propagatePositionsVelocities();
        void setBathVariables();
        double kineticEnergy;
        double kineticEnergyScaleFactor;
        double dt2;
        double dt4;
        double dt8;
    };
#endif

