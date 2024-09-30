#include "simpleModelDatabase.h"
/*! \file simpleModelDatabase.cpp */

/*!
Resize arrays and set whether to save various data entries
*/
simpleModelDatabase::simpleModelDatabase(int numberOfParticles, string fn, NcFile::FileMode mode,
                            bool saveVelocities, bool saveTypes, bool saveForces)
    :BaseDatabaseNetCDF(fn,mode)
    {
    velocity = saveVelocities;
    type = saveTypes;
    force = saveForces;

    N=numberOfParticles;
    dof = 3*N;
    val=0.0;
    vec.resize(dof);
    switch(mode)
        {
        case NcFile::read:
            GetDimVar();
            break;
        case NcFile::write:
            GetDimVar();
            break;
        case NcFile::replace:
            SetDimVar();
            break;
        case NcFile::newFile:
            SetDimVar();
            break;
        default:
            ;
        };
    };

void simpleModelDatabase::SetDimVar()
    {
    //Set the dimensions
    recDim = File.addDim("record");
    nDim = File.addDim("numberOfParticles", N);
    dofDim = File.addDim("spatialDegreesOfFreedom", dof);
    unitDim = File.addDim("unit",1);

    //Set the variables
    timeVar = File.addVar("time",ncDouble,recDim);
    positionVar = File.addVar("R3position",ncDouble,{recDim,dofDim});
    barycentricPositionVar = File.addVar("barycentricPosition",ncDouble,{recDim,dofDim});
    faceIndexVar = File.addVar("faceIndex",ncInt,{recDim,nDim});

    if(velocity)
        velocityVar = File.addVar("velocity",ncDouble,{recDim,dofDim});
    if(force)
        forceVar = File.addVar("force",ncDouble,{recDim,dofDim});
    if(type)
        typeVar = File.addVar("type",ncInt,{recDim,dofDim});
    }

void simpleModelDatabase::GetDimVar()
    {
    //Get the dimensions
    recDim = File.getDim("record");
    nDim = File.getDim("numberOfParticles");
    dofDim = File.getDim("spatialDegreesOfFreedom");
    unitDim = File.getDim("unit");

    //Get the variables
    timeVar = File.getVar("time");
    positionVar = File.getVar("R3position");
    barycentricPositionVar = File.getVar("barycentricPosition");
    faceIndexVar = File.getVar("faceIndex");

    if(velocity)
        velocityVar = File.getVar("velocity");
    if(force)
        forceVar = File.getVar("force");
    if(type)
        typeVar = File.getVar("type");
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
void simpleModelDatabase::writeState(STATE s, double time, int rec)
    {
    int record = rec;
    if(record<0)
        record = recDim.getSize();
    std::vector<double> r3pos(dof);
    std::vector<double> baryPos(dof);
    std::vector<double> vel(dof);
    std::vector<double> forceVector(dof);
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


    timeVar.putVar({record},&time);
    positionVar             .putVar({record,dofDim.getSize()},{1,dofDim.getSize()}, &r3pos[0]);
    barycentricPositionVar  .putVar({record,dofDim.getSize()},{1,dofDim.getSize()}, &baryPos[0]);
    faceIndexVar            .putVar({record,nDim.getSize()},{1,nDim.getSize()},&faceIdx[0]);
    if(velocity)
        velocityVar         .putVar({record,dofDim.getSize()},{1,dofDim.getSize()},&vel[0]);
    if(force)
        forceVar            .putVar({record,dofDim.getSize()},{1,dofDim.getSize()},&forceVector[0]);
    if(type)
        typeVar             .putVar({record,nDim.getSize()},{1,nDim.getSize()},&typeVector[0]);

    File.sync();
    };
/*!
This loads the targeted record into the state.  If you want to load the last
state in the file (or otherwise want to know what the options are), use the
GetNumRecs() function of the database.  This will return an int, and the final
record will be (GetNumRecs()-1).  Please also note that if you have not saved
velocity and force information in the database and then run some equation of
motion, the first timestep might be incorrect (i.e., with updaters that use a
partial-timestep-scheme that first updates a position via current velocity,
etc).
*/
void simpleModelDatabase::readState(STATE s, int rec)
    {
    int totalRecords = GetNumRecs();
    if (rec >= totalRecords)
        {
        printf("Trying to read a database entry that does not exist\n");
        throw std::exception();
        };

    std::vector<double> baryPos(dof,0.);
    std::vector<double> vel(dof,0.);
    std::vector<int> typeVector(N,0);
    std::vector<int> faceIdx(N,0);

    barycentricPositionVar.getVar({rec,dofDim.getSize()},{1,dofDim.getSize()},&baryPos[0]);

    faceIndexVar.getVar({rec,nDim.getSize()},{1,nDim.getSize()},&faceIdx[0]);
    if(velocity)
        velocityVar.getVar({rec,dofDim.getSize()},{1,dofDim.getSize()},&vel[0]);
    if(type)
        {
        s->types.resize(N);
        typeVar.getVar({rec,nDim.getSize()},{1,nDim.getSize()},&typeVector[0]);
        };


    for (int ii = 0; ii < N; ++ii)
        {
        point3 baryP(baryPos[3*ii],baryPos[3*ii+1],baryPos[3*ii+2]);
        vector3 v(vel[3*ii],vel[3*ii+1],vel[3*ii+2]);
        s->positions[ii] = meshPosition(baryP,faceIdx[ii]);
        if(velocity)
            s->velocities[ii] = v;
        if(type)
            s->types[ii] = typeVector[ii];
        };

    };

