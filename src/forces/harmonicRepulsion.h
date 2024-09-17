#ifndef harmonicRepulsion_H
#define harmonicRepulsion_H

#include "baseForce.h"

/*!
A force associated with a harmonic repulsive force between two points.  Given a
particle diameter \sigma, a stiffness k,
and a distance between points r, the energy is
E(r) = 0.5*k*(1-|r|/\sigma)^2
and the force is the gradient of this.  At the moment only monodisperse
interactions have been added, but changing this to read the radius of particles
is on the agenda
*/
class harmonicRepulsion : public force
    {
    public:
        harmonicRepulsion(double stiffness=1, double monodisperseRange=1, bool _monodisperse = true)
            {
            k=stiffness;
            sigma = monodisperseRange;
            monodisperse = _monodisperse;
            maximumInteractionRange = monodisperseRange;
            };

        virtual double pairwiseEnergy(vector3 separation,double distance);
        virtual vector3 pairwiseForce(vector3 separation,double distance);
    protected:
        double k=1;
        double sigma = 1;
        bool monodisperse = true;
    };
#endif
