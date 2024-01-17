#ifndef gaussianRepulsion_H
#define gaussianRepulsion_H

#include "baseForce.h"

class gaussianRepulsion : public force
    {
    public:
        gaussianRepulsion(double scale=1, double s=3)
            {
            alpha = scale;
            beta  = s;
            };

        virtual double pairwiseEnergy(vector3 separation,double distance);
        virtual vector3 pairwiseForce(vector3 separation,double distance);

    protected:
        double alpha;
        double beta;
    };
#endif
