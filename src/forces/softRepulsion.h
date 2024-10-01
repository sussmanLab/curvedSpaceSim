#ifndef softRepulsion_H
#define softRepulsion_H

#include "baseForce.h"

/*!
A force associated with a soft repulsive force between two points.  Given a
particle diameter \sigma, a stiffness k,
and a distance between points r, the energy is
E(r) = a^-1*k*(1-|r|/\sigma)^a
and the force is the gradient of this.  At the moment only monodisperse
interactions have been added, but changing this to read the radius of particles
is on the agenda
*/
class softRepulsion : public force
    {
    public:
        softRepulsion(double exponent = 2.5,double stiffness=1, double monodisperseRange=1, bool _monodisperse = true)
            {
            alpha = exponent;
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
        double alpha = 2.5; //default to Hertzian
        bool monodisperse = true;
    };
#endif
