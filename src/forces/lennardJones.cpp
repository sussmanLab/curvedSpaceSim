#include "lennardJones.h"

double lennardJones::potential(double distance)
    {
    return 4*eps*(pow(sigma/distance,12) - pow(sigma/distance, 6));
    }

double lennardJones::potentialDerivative(double distance)
    {
    return 4*eps*(6*pow(sigma/distance, 6)/distance - 12*pow(sigma/distance, 12)/distance); 
    }	    

double lennardJones::pairwiseEnergy(vector3 separation,double distance)
    {
    double ans = 0;
    if(monodisperse)
        {
	//implementing truncated + shifted lennard jones energy/force 
	//so that submesh cutoffs make sense
        if(distance < cutoffCoefficient*sigma)
	    {
	    ans = potential(distance) - potential(cutoffCoefficient*sigma);
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
        if(distance < cutoffCoefficient*sigma)
            {
            ans = (potentialDerivative(distance))*separation;
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
    model->findNeighbors(cutoffCoefficient*maximumInteractionRange);

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

double lennardJones::computeEnergy(bool verbose)
    {
    model->findNeighbors(cutoffCoefficient*maximumInteractionRange);
    energy = 0;
    for (int ii = 0; ii < model->N; ++ii)
        {
        int neighborNumber = model->neighbors[ii].size();
        for (int jj = 0; jj < neighborNumber; ++jj)
            energy += pairwiseEnergy(model->neighborVectors[ii][jj],model->neighborDistances[ii][jj]);
        }
    return energy;
    };
