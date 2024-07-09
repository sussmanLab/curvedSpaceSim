#ifndef DATABASENETCDF_H
#define DATABASENETCDF_H

#include <netcdfcpp.h>
#include "baseDatabase.h"

/*! \file DatabaseNetCDF.h */
//! A base class that implements a details-free  netCDF4-based data storage system. The unlimited dimension should always be named "record"
/*!
BaseDatabase just provides an interface to a file and a mode of operation.
*/
class BaseDatabaseNetCDF : public baseDatabase
    {
    public:
        //!The NcFile itself
        NcFile File;

        //!The default constructor starts a bland filename in readonly mode
        BaseDatabaseNetCDF(string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly);
        //!read the number of records in the database
        int GetNumRecs(){
                    NcDim *rd = File.get_dim("record");
                    return rd->size();
                    };
    };

#endif
