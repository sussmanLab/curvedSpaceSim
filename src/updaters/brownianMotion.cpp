#ifndef brownianMotion_H
#define brownianMotion_H

#include "baseUpdater.h"

// class like gradient descent, but with random kicks applied to the force 
// to give ``activity'' to particles

class brownianMotion : public updater
    {
    public:
        brownianMotion(double _dt)
            {
            deltaT = _dt;
            };

        virtual void performUpdate();

        //!A sample function, mostly to show how to use manipulateUpdaterData
        virtual double getMaxForce();
        //!A sample function, mostly to show how to use manipulateUpdaterData
        virtual double getForceNorm();

        double squaredTotalForceNorm;
        double maximumForceNorm;
    protected:
        vector<vector3> displacements;

    };
#endif
