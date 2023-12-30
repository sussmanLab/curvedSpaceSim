#ifndef gradientDescent_H
#define gradientDescent_H

#include "baseUpdater.h"


class gradientDescent : public updater
    {
    public:
        gradientDescent(double _dt)
            {
            deltaT = _dt;
            };  

        virtual void performUpdate();

        //!A sample function, mostly to show how to use manipulateUpdaterData
        virtual double getMaxForce();
        //!A sample function, mostly to show how to use manipulateUpdaterData
        virtual double getForceNorm();

    protected:
        vector<vector3> displacements;

    };
#endif
