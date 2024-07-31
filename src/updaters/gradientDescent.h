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

	//if we're switching between updaters (e.g. annealing), 
	//gradient descent needs to explicitly set that velocity transport is false 
	virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            model->particleShiftsRequireVelocityTransport = false;
	    initializeFromModel();
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
