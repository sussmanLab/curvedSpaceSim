#ifndef fireMinimization_H
#define fireMinimization_H

#include "velocityVerletNVE.h"

/*! \file fireMinimization.h */
//!Implement energy minimization via the FIRE algorithm
/*!
This class uses the "FIRE" algorithm to perform an energy minimization.
Each timestep, though, is a complete minimization (i.e., will run for the maximum number
of iterations, or until a target tolerance has been achieved, or whatever stopping condition the user
sets.
*/
class fireMinimization : public velocityVerletNVE
    {
    public:
        fireMinimization(){deltaT = 0.001;alpha = 0.99;};

        virtual void performUpdate(minimizeByFire());

        //! Minimize to the target force tolerance, or to the maximum number of iterations
        void minimizeByFire();

        void fireStep();

        //!set all of the fire parameters
        void setFIREParameters(int _maximumIterations,double _deltaT, double _alphaStart, double _deltaTMax, double _deltaTMin, double _deltaTInc, double _deltaTDec, double _alphaDec, int _nMin, double _forceCutoff, double _alphaMin = 0.75);

        double forceMax;
        double power;
        double forceNorm;
        double velocityNorm;
    protected:
        int maximumIterations = 1000;
        double alphaStart=0.99;
        double deltaTMax=0.1;
        double deltaTInc = 1.1;
        double deltaTMin = 1e-5;
        double deltaTDec = 0.95;
        double alphaDec = 0.9;
        double forceCutoff = 1e-12;
        double alphaMin = 0.0;
        double alpha;
        int nMin =4;
        int nSinceNegativePower=0;

        double updaterVectorDotProduct(vector<vector3> &v1,vector<vector3> &v2);
    };
#endif
