#ifndef velocityVerletNVE_H
#define velocityVerletNVE_H

#include "baseUpdater.h"

class velocityVerletNVE : public updater
    {
    public:
        velocityVerletNVE(double _dt)
            {
            deltaT = _dt;
            };
        virtual void performUpdate();

        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            model->particleShiftsRequireVelocityTransport = true;
            initializeFromModel();
            };
    protected:
        vector<vector3> displacements;
    };
#endif
