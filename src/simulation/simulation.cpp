#include "simulation.h"
/*! \file simulation.cpp */

/*!
Add a pointer to the vector of updaters, and give that updater a reference to the
model...
*/
void simulation::addUpdater(UpdaterPtr _upd, ConfigPtr _config)
    {
    _upd->setModel(_config);
    _upd->setSimulation(getPointer());
    updaters.push_back(_upd);
    };

/*!
Add a pointer to the vectors of force computers, and give that FC a reference to
the model...
*/
void simulation::addForce(ForcePtr _force, ConfigPtr _config)
    {
    _force->setModel(_config);
    forceComputers.push_back(_force);
    };

/*!
Set a pointer to the configuration (aka the model)
*/
void simulation::setConfiguration(ConfigPtr _config)
    {
    configuration = _config;
    };

/*!
pass a value of dt to all updaters associated with the simulation
*/
void simulation::setIntegrationTimestep(double dt)
    {
    integrationTimestep = dt;
    for (int u = 0; u < updaters.size(); ++u)
        {
        auto upd = updaters[u].lock();
        upd->setDeltaT(dt);
        };
    };

/*!
Loop over all associated updaters and set reproducible dynamics to the value of the boolean argument. For some updaters this will do nothing. For any updater that uses a random number generator, it will ask that updater to use a fixed (specific) seed.
*/
void simulation::setReproducible(bool reproducible)
    {
    for (int u = 0; u < updaters.size(); ++u)
        {
        auto upd = updaters[u].lock();
        upd->setReproducible(reproducible);
        };
    };

/*!
Loop over any force associated with the simulation, and call its computeForces
function.  For the first call, make sure we start with a zero force vector,
otherwise update the total force on each particle
*/
void simulation::computeForces()
    {
    auto Conf = configuration.lock();
    for (unsigned int f = 0; f < forceComputers.size(); ++f)
        {
        auto frc = forceComputers[f].lock();
        bool zeroForces = (f==0);
        frc->computeForces(Conf->forces,zeroForces);
        };
    };

/*!
Calls the configuration to displace the degrees of freedom
*/
void simulation::moveParticles(vector<vector3> &displacements)
    {
    auto Conf = configuration.lock();
    Conf->moveParticles(displacements);
    };

/*!
Call all updaters to advance the system one time step; for convenience, update the simulation variables to account for the timestep and time
*/
void simulation::performTimestep()
    {
    integerTimestep += 1;
    Time += integrationTimestep;

    //perform any updates, one of which should probably be an EOM
    for (int u = 0; u < updaters.size(); ++u)
        {
        auto upd = updaters[u].lock();
        upd->Update(integerTimestep);
        };
    };

/*!
Function to calculate the outer product, overloaded to work on both double[3] and vector3 objects.
Adds to an existing array of doubles (flattened, 9 entries). 
*/
void simulation::addOuterProduct(const double u[3], const double v[3], double (&result)[9]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[3*i+j] += u[i] * v[j];  
        }
    }
}

void simulation::addOuterProduct(const vector3 u, const vector3 v, double (&result)[9]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[3*i+j] += u[i] * v[j]; 
        }
    }
}

/* Compute the ``stress'' tensor of a monodisperse system using the virial approximation 
 * sigma = (1/2) rho kb <(vi outer vi)> + <f_ij outer dr_ij>/(2*d*V), d = 2 (for a surface). 
 * Modifies a (flattened) 3x3 array expressing the Euclidean stress tensor on the surface, "stress". 
 * Note that the tangent spaces of different particles pairs are different, and so some questions
 * about the interpretation of this global stress tensor remain. For the pressure (i.e., the 
 * tensor's trace) this is not an issue.
 */
void simulation::computeMonodisperseStress(vector<double> &stress) 
    {
    //note: locks are tools to prevent disparate threads from accessing the same
    //object simultaneously.  
    auto conf = configuration.lock();  
    double kb = 1; //currently always working in units of absolute temp, but altering kb here *does* change sim units
    auto f0 = forceComputers[0].lock(); 
    double particleSize = f0->maximumInteractionRange; //technically this uses the particle size from just the first force computer

    //calculate density, currently number/area so we don't break w/ multiple force computers & to match mass=1 theory
    double meshSurfaceArea = conf->space->getArea();  //space is subcomponent of model
    int Ndof = conf->N; 
    double density = Ndof/meshSurfaceArea; //using number/area density definition to avoid mismatch b/w differing cutoffs of force computers

    //positions and velocities will be necessary later
    vector<vector3> velocities = conf->velocities; 
    conf->fillEuclideanLocations(); 
    vector<double3> positions = conf->euclideanLocations;  
    //we need to calculate fij outer drij for *every* force computer
    conf->findNeighbors(particleSize); 
   
    //now, we need to calculate the second term for every force computer. This calculation takes a lot of 
    //its structure from computeForces within baseForce.cpp, and the necessary double loop will be pretty slow.
    //note that 9 and 3 below are 3*d and 1*d respectively, where d is the dimension of the force/sep/velocity vectors. they're euclidean!  
    double fOuterDR [9] =  {0,0,0,0,0,0,0,0,0};      
    double vOuterv [9] = {0,0,0,0,0,0,0,0,0};
    bool calcvOuterv = true; //only want to calc v outer v during first force computer, then never again (tm)
    
    for (unsigned int f = 0; f < forceComputers.size(); ++f)
        {
        auto frc = forceComputers[f].lock();
	
	      for (int ii = 0; ii < Ndof; ++ii)
            { 
	          int neighborNumber = conf->neighbors[ii].size();
	          if (calcvOuterv) addOuterProduct(velocities[ii], velocities[ii], vOuterv); //placed here for loop convenience, only calculate during first force computer;
	          
            for (int jj = 0; jj < neighborNumber; ++jj)
	              {
                vector3 separation = conf->neighborVectors[ii][jj]; //specifies the separation vector between particle i and its jth neighbor
                double distance = conf->neighborDistances[ii][jj]; //specifies the distance between particle i and its jth neighbor
                vector3 separationWithDistance = distance*separation;
                vector3 force = -frc->pairwiseForce(separation, distance); //calculates the associated force
                //note: negative above is a consequence of the sign convention we choose for force and separation. 
                //The standard formula uses positive force. 
	              //now add each outer product between force and separation in turn into fOuterDR
                addOuterProduct(force, separationWithDistance, fOuterDR);     
	              }
            }
            
        calcvOuterv = false; //after first force computer loop, we're done calculating v outer v 	
        }

    //now that we've actually calculated v(outer)v and f(outer)dr, we can use them to calculate stress
    //sigma = (1/2) rho <(vi outer vi)> + (f_ij outer dr_ij)/(2*d*V), d = 2 -- dimension here is NOT the same as above, as now we are restricting to surface
    for (int ind = 0; ind < 9; ind ++) 
        {
        vOuterv[ind] = density*vOuterv[ind]/(2*Ndof); //trace of vouterv/(2*Ndof) is average KE (aka temperature) 
        fOuterDR[ind] = fOuterDR[ind]/(4*meshSurfaceArea); //(1/4) here is (1/2) for double counts, (1/2) for dimension
        stress[ind] = vOuterv[ind]+fOuterDR[ind]; 
        }
    }
