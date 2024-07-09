#include "simpleModelDatabase.h"
/*! \file simpleModelDatabase.cpp */

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
        case NcFile::ReadOnly:
            GetDimVar();
            break;
        case NcFile::Write:
            GetDimVar();
            break;
        case NcFile::Replace:
            SetDimVar();
            break;
        case NcFile::New:
            SetDimVar();
            break;
        default:
            ;
        };
    };

void simpleModelDatabase::SetDimVar()
    {
    //Set the dimensions
    recDim = File.add_dim("record");
    nDim = File.add_dim("numberOfParticles", N);
    dofDim = File.add_dim("spatialDegreesOfFreedom", dof);
    unitDim = File.add_dim("unit",1);

    //Set the variables
    timeVar = File.add_var("time",ncDouble,recDim);
    positionVar = File.add_var("R3position",ncDouble,recDim,dofDim);
    barycentricPositionVar = File.add_var("barycentricPosition",ncDouble,recDim,dofDim);
    faceIndexVar = File.add_var("faceIndex",ncInt,recDim,nDim);

    if(velocity)
        velocityVar = File.add_var("velocity",ncDouble,recDim,dofDim);
    if(force)
        forceVar = File.add_var("force",ncDouble,recDim,dofDim);
    if(type)
        typeVar = File.add_var("type",ncInt,recDim,dofDim);
    }

void simpleModelDatabase::GetDimVar()
    {
    //Get the dimensions
    recDim = File.get_dim("record");
    nDim = File.get_dim("numberOfParticles");
    dofDim = File.get_dim("spatialDegreesOfFreedom");
    unitDim = File.get_dim("unit");

    //Get the variables
    timeVar = File.get_var("time");
    positionVar = File.get_var("R3position");
    barycentricPositionVar = File.get_var("barycentricPosition");
    faceIndexVar = File.get_var("faceIndex");

    if(velocity)
        velocityVar = File.get_var("velocity");
    if(force)
        forceVar = File.get_var("force");
    if(type)
        typeVar = File.get_var("type");
    }

void simpleModelDatabase::writeState(STATE s, double time, int rec)
    {
    int record = rec;
    if(record<0)
        record = recDim->size();
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
        }


    timeVar  -> put_rec(&time, record);
    positionVar             ->put_rec(&r3pos[0],record);
    barycentricPositionVar  ->put_rec(&baryPos[0],record);
    faceIndexVar            ->put_rec(&faceIdx[0],record);
    if(velocity)
        velocityVar         ->put_rec(&vel[0],record);
    if(force)
        forceVar            ->put_rec(&forceVector[0],record);
    if(type)
        typeVar             ->put_rec(&typeVector[0],record);

    File.sync();
    };

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

    barycentricPositionVar->set_cur(rec);
    barycentricPositionVar -> get(&baryPos[0],1,dofDim->size());

    faceIndexVar->set_cur(rec);
    faceIndexVar-> get(&faceIdx[0],1,nDim->size());
    if(velocity)
        {
        velocityVar ->set_cur(rec);
        velocityVar ->get(&vel[0],1,dofDim->size());
        };
    if(type)
        {
        s->types.resize(N);
        typeVar ->set_cur(rec);
        typeVar-> get(&typeVector[0],1,nDim->size());
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

