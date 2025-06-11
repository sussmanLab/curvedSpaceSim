#include "nonInteractingForce.h"

double nonInteractingForce::pairwiseEnergy(vector3 separation, double distance) 
    {
    return zeroPoint;
    } 

vector3 nonInteractingForce::pairwiseForce(vector3 separation, double distance) 
    {
    return vector3(0,0,0);
    } 
