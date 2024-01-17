#include "rieszRepulsion.h"

double rieszRepulsion::pairwiseEnergy(vector3 separation,double distance)
    {
    return alpha/pow(distance,beta);
    };

vector3 gaussianRepulsion::pairwiseForce(vector3 separation,double distance)
    {
    double prefactor = beta*alpha/pow(distance,beta+1);
    return -prefactor*separation;
    };
