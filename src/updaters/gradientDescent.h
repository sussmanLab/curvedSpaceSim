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

        virtual void performUpdater();

        virtual double getMaxForce();

    protected:
        vector<vector3> displacements;

    }
