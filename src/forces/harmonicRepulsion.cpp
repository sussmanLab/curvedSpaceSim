#include "harmonicRepulsion.h"

double harmonicRepulsion::pairwiseEnergy(vector3 separation,double distance)
    {
    double ans = 0;
    if(monodisperse)
        {
        if(distance < sigma)
            ans = 0.5*k*(1-distance/sigma)*(1-distance/sigma);
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
            //kDelta is dU/d\delta
            double kDelta = k*(1.0-distance/sigma);
            ans = -(kDelta/sigma)*separation;
            }
        }
    else
        {
        UNWRITTENCODE("non-monodisperse harmonic repulsions not implemented");
        }
    return ans;
    };
