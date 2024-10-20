#include "lennardJones.h"

double lennardJones::pairwiseEnergy(vector3 separation,double distance)
    {
    double ans = 0;
    if(monodisperse)
        {
	//implementing truncated + shifted lennard jones energy/force 
	//so that submesh cutoffs make sense
        if(distance < 2*sigma)
	    {
            //subtracting u_lj(r=2*sigma) here
	    double uOfTwoSigma = -63.0/4096.0; 
	    ans = 4*eps*((pow(sigma/distance,12) - pow(sigma/distance, 6))  - uOfTwoSigma );
            }
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
	    double FofTwoSigma = -(93.0/2048.0)*(1/sigma);
            ans = 4*eps*(6*pow(sigma/distance,6)/distance - 12*pow(sigma/distance,12)/distance - FofTwoSigma)*separation;
            }
        }
    else
        {
        UNWRITTENCODE("non-monodisperse l-j interactions not implemented");
        }
    return ans;
    };

//special compute forces so we double the interaction range relative to particle size
void lennardJones::computeForces(vector<vector3> &forces,bool zeroOutForce, int type)
    {
    if(forces.size() != model->N)
        forces.resize(model->N);
    model->findNeighbors(2*maximumInteractionRange);

    for (int ii = 0; ii < model->N; ++ii)
        {
        if(zeroOutForce)
            forces[ii] = vector3(0.0,0.0,0.0);
        int neighborNumber = model->neighbors[ii].size();
        for (int jj = 0; jj < neighborNumber; ++jj)
            {
            forces[ii] += pairwiseForce(model->neighborVectors[ii][jj],model->neighborDistances[ii][jj]);
            }
        }
    };
