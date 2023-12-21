#include "simulation.h"
/*! \file simulation.cpp */

/*!
Add a pointer to the list of updaters, and give that updater a reference to the
model...
*/
void simulation::addUpdater(UpdaterPtr _upd, ConfigPtr _config)
    {
    _upd->setModel(_config);
    _upd->setSimulation(getPointer());
    updaters.push_back(_upd);
    };

/*!
Add a pointer to the list of force computers, and give that FC a reference to the
model...
*/
void simulation::addForce(ForcePtr _force, ConfigPtr _config)
    {
    _force->setModel(_config);
    forceComputers.push_back(_force);
    };

/*!
Set a pointer to the configuration
*/
void simulation::setConfiguration(ConfigPtr _config)
    {
    configuration = _config;
    };

/*!
\post the cell configuration and e.o.m. timestep is set to the input value
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
\post the updaters are set to be reproducible if the boolean is true, otherwise the RNG is initialized
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
Calls all force computers, and evaluate the self force calculation if the model demands it
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
Call all relevant functions to advance the system one time step; every sortPeriod also call the
spatial sorting routine.
\post The simulation is advanced one time step
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
