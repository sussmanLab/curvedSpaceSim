#ifndef gradientDescent_H
#define gradientDescent_H

#include "baseUpdater.h"
//! A class which performs a bog-standard gradient descent 
class gradientDescent : public updater
    {
    public:
        gradientDescent(double _dt)
            {
            deltaT = _dt;
            };  

        virtual void performUpdate();
    protected:
        vector<vector3> displacements;

    };
#endif
