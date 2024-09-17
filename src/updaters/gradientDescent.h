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

	//if we're switching between updaters (e.g. annealing), 
	//gradient descent needs to explicitly set that velocity transport is false 
	virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
	    initializeFromModel();
            };

        virtual void performUpdate();
    protected:
        vector<vector3> displacements;

    };
#endif
