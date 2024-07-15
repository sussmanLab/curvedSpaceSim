#include "lennardJones.h"

double lennardJones::pairwiseEnergy(vector3 separation,double distance)
    {
    double ans = 0;
    if(monodisperse)
        {
	//implementing truncated + shifted lennard jones energy/force 
	//so that submesh cutoffs make sense
        if(distance < 2*sigma)
            //subtracting u_lj(r=2*sigma) here
	    double uof2sigma = -63.0/4096.0; 
	    ans = 4*eps*((pow(sigma/distance,12) - pow(sigma/distance, 6))  - uof2sigma );
        }
    else
        {
        UNWRITTENCODE("non-monodisperse l-j interactions not implemented");
        }
    return ans;
    };

//note that computeForces expects neighbor vectors to be *unit vectors*   
vector3 lennardJones::pairwiseForce(vector3 separation,double distance)
    {
    vector3 ans = vector3(0.0,0.0,0.0);
    if(monodisperse)
        {
        if(distance <= 2*sigma)
            {
            //f = -grad(U), but the negative is in the updater rather than here
	    //also -- still truncating and shifting force so that it's zero at the cutoff
	    double Fof2sigma = -(93.0/2048.0)*(1/sigma);
            ans = 4*eps*(6*pow(sigma/distance,6)/distance - 12*pow(sigma/distance,12)/distance - Fof2sigma)*separation;
            }
        }
    else
        {
        UNWRITTENCODE("non-monodisperse l-j interactions not implemented");
        }
    return ans;
    };
