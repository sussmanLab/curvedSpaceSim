#include "gaussianRepulsion.h"

double gaussianRepulsion::pairwiseEnergy(vector3 separation,double distance)
    {
    return alpha*exp(-distance*distance/twoSigmaSquared)/(sigmaSqrtTwoPi);
    };

vector3 gaussianRepulsion::pairwiseForce(vector3 separation,double distance)
    {
    double prefactor = distance*alpha* exp(-distance*distance/twoSigmaSquared) / (sigmaThreeHalvesSqrtTwoPi);
    return -prefactor*separation;
    };
