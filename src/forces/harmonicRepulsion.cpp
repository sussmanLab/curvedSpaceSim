#include "harmonicRepulsion.h"

double harmonicRepulsion::pairwiseEnergy(vector3 separation,double distance)
    {
    double ans = 0;
    if(monodisperse)
        {
        if(distance < sigma)
            ans = 0.5*k*(sigma-distance)*(sigma-distance);
        }
    else
        {
        UNWRITTENCODE("non-monodisperse harmonic repulsions not implemented");
        }
    return ans;
    };

//note that computeForces expects neighbor vectors to be *unit vectors*   
vector3 harmonicRepulsion::pairwiseForce(vector3 separation,double distance)
    {
    vector3 ans = vector3(0.0,0.0,0.0);
    if(monodisperse)
        {
        if(distance <= sigma)
            {
            ans = -k*(sigma-distance)*separation;
            }
        }
    else
        {
        UNWRITTENCODE("non-monodisperse harmonic repulsions not implemented");
        }
    return ans;
    };
