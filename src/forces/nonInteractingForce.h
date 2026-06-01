#ifndef nonInteractingForce_H
#define nonInteractingForce_H

#include "baseForce.h"

/*!
A force corresponding to constant energy, intended for noninteracting particles.
*/
class nonInteractingForce : public force
    {
    public:
	nonInteractingForce(double basicValue = 0)
            {
	    zeroPoint = basicValue;
	    }
	virtual double pairwiseEnergy(vector3 separation, double distance);
        virtual vector3 pairwiseForce(vector3 separation, double distance);

    protected:
       double zeroPoint=0;	
    };
#endif
