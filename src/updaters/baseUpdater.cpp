#include "baseUpdater.h"
/*! \file baseUpdater.cpp" */

void updater::performUpdate()
    {
    cout << "in the base updater... that's odd..." << endl;
    sim->computeForces();
    };

int updater::getNTotal()
    {
    if(updaterData.size()<1)
        updaterData.resize(1);
    nTotal = model ->getNumberOfParticles();

    return nTotal;
    }
