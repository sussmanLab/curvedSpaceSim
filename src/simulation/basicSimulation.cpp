#include "basicSimulation.h"
/*! \file basicSimulation.cpp */

/*!
Initialize all of the shared points, set default values of things
*/
basicSimulation::basicSimulation(): integerTimestep(0), Time(0.),integrationTimestep(0.01)
    {
    };

/*!
\post the cell configuration and e.o.m. timestep is set to the input value
*/
void basicSimulation::setCurrentTime(double _cTime)
    {
    Time = _cTime;
    //auto Conf = cellConfiguration.lock();
    //Conf->setTime(Time);
    };


