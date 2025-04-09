#include "simpleModelDatabase.h"
#include "baseDatabase.h"
#include "baseHDF5Database.h"
#include "debuggingHelp.h"
/*! \file simpleModelDatabase.cpp */

/*!
Resize arrays and set whether to save various data entries
*/
simpleModelDatabase::simpleModelDatabase(int numberOfParticles, string fn, fileMode::Enum _mode,
                            bool saveVelocities, bool saveTypes, bool saveForces)
    :baseHDF5Database(fn,mode)
    {
    velocity = saveVelocities;
    type = saveTypes;
    force = saveForces;

    N = numberOfParticles;
    timeVector.resize(1);
    intVector.resize(N);
    doubleVector.resize(N);
    coordinateVector.resize(3*N);

    if(mode == fileMode::replace)
        {
        registerDatasets();
        };
    if(mode == fileMode::readwrite)
        {
        if (currentNumberOfRecords() ==0)
            registerDatasets();
        };
    logMessage(logger::verbose, "modelDatabase initialized");
    }

void simpleModelDatabase::registerDatasets()
    {
    registerExtendableDataset<double>("time", 1);
    registerExtendableDataset<double>("R3position", 3*N);
    registerExtendableDataset<double>("barycentricPosition", 3*N);
    registerExtendableDataset<int>("faceIndex", N);

    if(velocity)
        registerExtendableDataset<double>("velocity", 3*N);
    if(force)
        registerExtendableDataset<double>("force", 3*N);
    if(type)
        registerExtendableDataset<int>("type", N);
    }


unsigned long simpleModelDatabase::currentNumberOfRecords()
    {
    return getDatasetDimensions("time");
    }


/*!
Save a record that always saves particle positions and optionally a few other
data structures to the file Information is automatically extracted from the
member variables in the STATE.  Additionally, a scalar value will be saved with
the record (typically the current simulation time, but could be anything).  By
default (rec = -1) the record will be appended as a new record in the netcdf
file.  If a specific index is given, it will be saved to that record instead
(potentially overwriting existing information)
*/
void simpleModelDatabase::writeState(STATE s, double time, int record)
    {
    if(record >= 0)
        ERRORERROR("overwriting specific records not implemented at the moment");

    timeVector[0] = time;
    extendDataset("time", timeVector);


    std::vector<double> r3pos(3*N);
    std::vector<double> baryPos(3*N);
    std::vector<double> vel(3*N);
    std::vector<double> forceVector(3*N);

    std::vector<int> typeVector(N);
    std::vector<int> faceIdx(N);

    //be able to save both barycentric data (for loading) and euclidean particle positions (for convenience)
    s->fillEuclideanLocations();

    for(int ii = 0; ii < N; ++ii)
        {
        r3pos[3*ii+0] = s->euclideanLocations[ii].x;
        r3pos[3*ii+1] = s->euclideanLocations[ii].y;
        r3pos[3*ii+2] = s->euclideanLocations[ii].z;
        
        baryPos[3*ii+0] = s->positions[ii].x[0];
        baryPos[3*ii+1] = s->positions[ii].x[1];
        baryPos[3*ii+2] = s->positions[ii].x[2];
        
        faceIdx[ii] = s->positions[ii].faceIndex;

        if(velocity)
            {
            vel[3*ii+0] = s->velocities[ii][0];
            vel[3*ii+1] = s->velocities[ii][1];
            vel[3*ii+2] = s->velocities[ii][2];
            }
        if(force)
            {
            forceVector[3*ii+0] = s->forces[ii][0];
            forceVector[3*ii+1] = s->forces[ii][1];
            forceVector[3*ii+2] = s->forces[ii][2];
            }
        if(type)
            typeVector[ii] = s->types[ii];
        };

    extendDataset("R3position", r3pos);
    extendDataset("barycentricPosition", baryPos);
    extendDataset("faceIndex", faceIdx);

    if(velocity)
        extendDataset("velocity", vel);
    if(force)
        extendDataset("force", forceVector);
    if(type)
        extendDataset("type", typeVector);
    }

void simpleModelDatabase::readState(STATE s, int rec)
    {
    readDataset("time",timeVector,rec);


    std::vector<double> baryPos(3*N,0.);
    std::vector<double> vel(3*N,0.);
    std::vector<double> forceVector(3*N,0.);
    std::vector<int> typeVector(N,0);
    std::vector<int> faceIdx(N,0);

    readDataset("barycentricPosition",baryPos);
    readDataset("faceIndex",faceIdx);
    if(velocity)
        readDataset("velocity",vel);
    if(force)
        readDataset("force",forceVector);
    if(type)
        readDataset("type",typeVector);

    for (int ii = 0; ii < N; ++ii)
        {
        point3 baryP(baryPos[3*ii],baryPos[3*ii+1],baryPos[3*ii+2]);
        s->positions[ii] = meshPosition(baryP,faceIdx[ii]);
        if(velocity)
            {
            vector3 v(vel[3*ii],vel[3*ii+1],vel[3*ii+2]);
            s->velocities[ii] = v;
            }
        if(force)
            {
            vector3 f(forceVector[3*ii],forceVector[3*ii+1],forceVector[3*ii+2]);
            s->forces[ii] = f;
            }
        if(type)
            s->types[ii] = typeVector[ii];
        };

    };

