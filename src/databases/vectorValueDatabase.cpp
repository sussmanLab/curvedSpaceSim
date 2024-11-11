#include "vectorValueDatabase.h"
/*! \file vectorValueDatabase.cpp */

vectorValueDatabase::vectorValueDatabase(int vectorLength, string fn, NcFile::FileMode mode)
    :BaseDatabaseNetCDF(fn,mode)
    {
    N=vectorLength;
    val=0.0;
    vec.resize(N);
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

void vectorValueDatabase::SetDimVar()
    {
    //Set the dimensions
    recDim = File.addDim("record");
    dofDim = File.addDim("dof", N);
    unitDim = File.addDim("unit",1);

    //Set the variables
    vecVar = File.addVar("vector", ncDouble,{recDim, dofDim});
    valVar = File.addVar("value",ncDouble,recDim);
    }

void vectorValueDatabase::GetDimVar()
    {
    //Get the dimensions
    recDim = File.getDim("record");
    dofDim = File.getDim("dof");
    unitDim = File.getDim("unit");
    //Get the variables
    vecVar = File.getVar("vector");
    valVar = File.getVar("value");
    }

void vectorValueDatabase::writeState(vector<double> &_vec,double _val)
    {
    int rec = recDim.getSize();
    valVar.putVar({rec},&_val);
    vecVar.putVar({rec,0},{1,dofDim.getSize()},&_vec[0]);
    File.sync();
    };

void vectorValueDatabase::readState(int rec)
    {
    int totalRecords = GetNumRecs();
    if (rec >= totalRecords)
        {
        printf("Trying to read a database entry that does not exist\n");
        throw std::exception();
        };
        vecVar.getVar({rec,N},{1,dofDim.getSize()},&vec[0]);
        valVar.getVar({rec},{1}, &val);
    };
