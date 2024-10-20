#ifndef lennardJones_H
#define lennardJones_H

#include "baseForce.h"

class lennardJones : public force
    {
    public:
        lennardJones(double energyScale=1, double monodisperseRange=1, bool _monodisperse = true)
            {
	    cout << "creating lennard jones, energy scale " << energyScale << ", sigma " << monodisperseRange << endl;
            eps = energyScale;
            monodisperse = _monodisperse;
            sigma = monodisperseRange;
            };

        virtual double pairwiseEnergy(vector3 separation,double distance);
        virtual vector3 pairwiseForce(vector3 separation,double distance);
	virtual void computeForces(vector<vector3> &forces, bool zeroOutForce, int type);
    protected:
        double eps=1;
        double sigma = 1;
        bool monodisperse = true;
    };
#endif

