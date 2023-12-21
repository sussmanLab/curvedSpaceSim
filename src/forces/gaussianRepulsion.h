#ifndef gaussianRepulsion_H
#define gaussianRepulsion_H

#include "baseForce.h"

class gaussianRepulsion : public force
    {
    public:
        gaussianRepulsion(double strength, double variance)
            {
            alpha =strength;
            sigma = variance;
            twoSigmaSquared = 2.0*sigma*sigma;
            sqrtSigma = sqrt(sigma);
            sigmaSqrtTwoPi = sqrtTwoPi* sigma;
            sigmaThreeHalvesSqrtTwoPi=sigmaSqrtTwoPi*sqrtSigma;
            };

        virtual double pairwiseEnergy(vector3 separation,double distance);
        virtual vector3 pairwiseForce(vector3 separation,double distance);

    protected:
        double alpha;
        double sigma;
        double sqrtTwoPi = 2.50662827463100050241576528481104525300698674061;
        double twoSigmaSquared;
        double sigmaSqrtTwoPi;
        double sqrtSigma;
        double sigmaThreeHalvesSqrtTwoPi;
    };
#endif
