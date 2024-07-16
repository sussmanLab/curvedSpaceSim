#include "noseHooverNVT.h"

noseHooverNVT::noseHooverNVT(double _dt, double _T, double _tau, int _M)
    {
    deltaT = _dt;
    dt2 = 0.5*deltaT;
    dt4 = 0.25*deltaT;
    dt8 = 0.125*deltaT;
    chainLength = _M;
    temperature = _T;
    tau = _tau;
    bathVariables.resize(chainLength+1);

    double4 zeroVector; zeroVector.x =0.0; zeroVector.y=0.0; zeroVector.z=0.0; zeroVector.w=0.0;
    for (int ii = 0; ii < bathVariables.size(); ++ii)
        bathVariables[ii] = zeroVector;
    };

/*!
Additionally, use the observation in the Mol Phys paper to set the masses of the chain of thermostats
*/
void noseHooverNVT::setBathVariables()
    {
    bathVariables[0].w = 2.0*(Ndof-2)*temperature*tau*tau;
    for (int ii = 1; ii < bathVariables.size(); ++ii)
        bathVariables[ii].w = temperature*tau*tau;

    kineticEnergy = bathVariables[0].w;
    kineticEnergyScaleFactor = 1.0;
    }

/*!
The implementation here closely follows algorithms 30 - 32 in Frenkel & Smit, generalized to the
case where the chain length is not necessarily always 2
*/
void noseHooverNVT::performUpdate()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
    
    //propagate the chain and rescale particle velocities
    propagateChain();
    for (int ii = 0; ii < Ndof; ++ii)
        model->velocities[ii] = kineticEnergyScaleFactor*model->velocities[ii];

    //propagate the positions and velocities of particles
    propagatePositionsVelocities();

    //propagate the chain and rescale particle velocities
    propagateChain();
    for (int ii = 0; ii < Ndof; ++ii)
        model->velocities[ii] = kineticEnergyScaleFactor*model->velocities[ii];
    };


/*!
This part of the algorithm partially updates the chain positions and velocities. It should be
called twice per time step
*/
void noseHooverNVT::propagateChain()
    {
    double exponentialFactor=0;
    //partially update bath velocities and accelerations over a quarter-timestep
    for(int ii = chainLength-1;ii>0;--ii)
        {
        //update the acceleration: G = (Q_{i-1}*v_{i-1}^2 - T)/Q_i
        bathVariables[ii].z = (bathVariables[ii-1].w*bathVariables[ii-1].y*bathVariables[ii-1].y-temperature)/bathVariables[ii].w;
        //the exponential factor is exp(-dt*v_{i+1}/2)
        exponentialFactor = exp(-dt8*bathVariables[ii+1].y);
        bathVariables[ii].y *= exponentialFactor;
        bathVariables[ii].y += bathVariables[ii].z*dt4;
        bathVariables[ii].y *= exponentialFactor;
        }
    bathVariables[0].z = (2.0*kineticEnergy - 2.0*(Ndof-2)*temperature)/bathVariables[0].w;
    exponentialFactor = exp(-dt8*bathVariables[1].y);
    bathVariables[0].y *= exponentialFactor;
    bathVariables[0].y += bathVariables[0].z*dt4;
    bathVariables[0].y *= exponentialFactor;

    //update the bath positions over a half timestep
    for (int ii = 0; ii < chainLength; ++ii)
        bathVariables[ii].x += dt2*bathVariables[ii].y;

    //get the scale factor for particle velocities and compute the KE
    kineticEnergyScaleFactor = exp(-dt2*bathVariables[0].y);
    kineticEnergy = kineticEnergyScaleFactor*kineticEnergyScaleFactor*kineticEnergy;

    //update the other quarter-timestep for the bath velocities and accelerations
    bathVariables[0].z = (2.0*kineticEnergy - 2.0*(Ndof-2)*temperature)/bathVariables[0].w;
    exponentialFactor = exp(-dt8*bathVariables[1].y);
    bathVariables[0].y *= exponentialFactor;
    bathVariables[0].y += bathVariables[0].z*dt4;
    bathVariables[0].y *= exponentialFactor;
    for (int ii = 1; ii < chainLength; ++ii)
        {
        bathVariables[ii].z = (bathVariables[ii-1].w*bathVariables[ii-1].y*bathVariables[ii-1].y-temperature)/bathVariables[ii].w;
        //the exponential factor is exp(-dt*v_{i+1}/2)
        exponentialFactor = exp(-dt8*bathVariables[ii+1].y);
        bathVariables[ii].y *= exponentialFactor;
        bathVariables[ii].y += bathVariables[ii].z*dt4;
        bathVariables[ii].y *= exponentialFactor;
        };
    }

/*!
This part of the algorithm updates positions and velocities of particles,
and requires a force calculation
*/
void noseHooverNVT::propagatePositionsVelocities()
    {
    kineticEnergy = 0.0;

    //move particles for the first half of the timestep, and compute the forces at the half timestep
    
    for (int ii = 0; ii < Ndof; ++ii)
        {
        displacements[ii] = dt2*model->velocities[ii];
        }
    sim->moveParticles(displacements);
    sim->computeForces();

    //update the positions and velocities (and evaluate the total KE) over the second half of the timestep

    for  (int ii = 0; ii < Ndof; ++ii)
        {
        model->velocities[ii] = model->velocities[ii] + (deltaT/model->masses[ii])*model->forces[ii];
        displacements[ii] = dt2*model->velocities[ii];
        double vDotV = model->velocities[ii]*model->velocities[ii];
        kineticEnergy += 0.5*(model->masses[ii])*vDotV;
        }
    sim->moveParticles(displacements);
    }


double noseHooverNVT::getTemperatureFromKE()
    {
    std::vector<double> kineticEnergies; 
    kineticEnergies.reserve(Ndof); 
    for (int i = 0; i < Ndof; i++) 
        {
	vector3 vel = model->velocities[i]; 
	kineticEnergies.push_back(vel*vel);
        }
    return (1/(2*Ndof))*std::accumulate(kineticEnergies.begin(), kineticEnergies.end(),0);

    }

