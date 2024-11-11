#include "softRepulsion.h"

double softRepulsion::pairwiseEnergy(vector3 separation,double distance)
    {
    double ans = 0;
    if(monodisperse)
        {
        if(distance < sigma)
            ans = (k/alpha)*pow((1-distance/sigma),alpha);
        }
    else
        {
        UNWRITTENCODE("non-monodisperse soft repulsions not implemented");
        }
    return ans;
    };

//note that computeForces expects neighbor vectors to be *unit vectors*   
vector3 softRepulsion::pairwiseForce(vector3 separation,double distance)
    {
    vector3 ans = vector3(0.0,0.0,0.0);
    if(monodisperse)
        {
        if(distance <= sigma)
            {
            //kDelta is dU/d\delta
            double kDelta = k*pow((1.0-distance/sigma),alpha-1);
            ans = -(kDelta/sigma)*separation;
            }
        }
    else
        {
        UNWRITTENCODE("non-monodisperse soft repulsions not implemented");
        }
    return ans;
    };
