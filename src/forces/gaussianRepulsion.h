#ifndef gaussianRepulsion_H
#define gaussianRepulsion_H

#include "baseForce.h"

/*!
A force associated with a gaussian repulsive force between two points.  Given a
gaussian of variance sigma, some prefactor alpha, and a vector distance between
points r, the energy is
E(r) = \frac{\alpha}{\sqrt{2\pi\sigma^2}} \exp{\frac{-r^2}{2\sigma^2}},
and the force is the gradient of this
*/
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
        //these functions are called all the time, so we precompute a few quantities
        double sqrtTwoPi = 2.50662827463100050241576528481104525300698674061;
        double twoSigmaSquared;
        double sigmaSqrtTwoPi;
        double sqrtSigma;
        double sigmaThreeHalvesSqrtTwoPi;
    };
#endif
