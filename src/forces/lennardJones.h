#ifndef lennardJones_H
#define lennardJones_H

#include "baseForce.h"

class lennardJones : public force
    {
    public:
        lennardJones(double energyScale=1, double monodisperseRange=1, bool _monodisperse = true)
            {
            eps = energyScale;
            monodisperse = _monodisperse;
            maximumInteractionRange = monodisperseRange;
            };

        virtual double pairwiseEnergy(vector3 separation,double distance);
        virtual vector3 pairwiseForce(vector3 separation,double distance);
    protected:
        double eps=1;
        double sigma = 1;
        bool monodisperse = true;
    };
#endif
