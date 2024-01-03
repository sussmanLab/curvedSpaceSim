#ifndef harmonicRepulsion_H
#define harmonicRepulsion_H

#include "baseForce.h"

class harmonicRepulsion : public force
    {
    public:
        harmonicRepulsion(double stiffnesss=1, double monodisperseRange=1, bool _monodisperse = true)
            {
            k=stiffness;
            sigma = monodisperseRange;
            monodisperse = _monodisperse;
            };

        virtual double pairwiseEnergy(vector3 separation,double distance);
        virtual vector3 pairwiseForce(vector3 separation,double distance);
    protected:
        double k=1;
        double sigma = 1;
        bool monodisperse = true;
    };
#endif
