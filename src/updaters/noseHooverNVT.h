#ifndef noseHooverNVT_H
#define noseHooverNVT_H

#include "baseUpdater.h"

/*! \file noseHooverNVT.h */
//! Implements NVT dynamics according to the Nose-Hoover equations of motion with a chain of thermostats
/*!
 This allows one to do standard NVT simulations.  A chain (whose length can be
 specified by the user) of thermostats is used to maintain the target
 temperature.  We closely follow the Frenkel & Smit update scheme, which is
 itself based on:
 Martyna, Tuckerman, Tobias, and Klein, Mol. Phys. 87, 1117 (1996) */

class noseHooverNVT : public updater
    {
    public:
        //here tau is the inverse of the frequency of the chain degrees of freedom, and M is the length of the NH chain
        noseHooverNVT(double _dt, double _T, double _tau=1.0, int _M=2);

        //!perform the MTTK symplectic NVT update
        virtual void performUpdate();

        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            model->particleShiftsRequireVelocityTransport = true;
            initializeFromModel();
            setBathVariables();
            };

        //!compute the instantaneous kinetic temperature
        double getTemperatureFromKE();

    protected:
        vector<vector3> displacements;
        double temperature;
        int chainLength;
        //!bathVariables contains (pos, vel, accel, mass) of the chain variables
        vector<double4> bathVariables;
        double tau;
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

