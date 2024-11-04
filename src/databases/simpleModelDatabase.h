#ifndef simpleModelDatabase_h
#define simpleModelDatabase_h

#include "DatabaseNetCDF.h"
#include "simpleModel.h"

/*! \file simpleModelDatabase.h */
//! A database to save states of simpleModels
/*!
There is one unlimited dimension, and each record can store a time and vector information associated with the model state.
The database always saves position information; flags on the constructor control whether other information is stored
*/
class simpleModelDatabase : public BaseDatabaseNetCDF
    {
    public:
        simpleModelDatabase(int numberOfParticles, string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly,
                            bool saveVelocities = true, bool saveTypes = false, bool saveForces = true);
        ~simpleModelDatabase(){File.close();};

        //! NcDims we'll use
        NcDim *recDim, *nDim,*dofDim, *unitDim;
        //! NcVars
        NcVar *timeVar, *positionVar, *barycentricPositionVar,*faceIndexVar, *velocityVar, *typeVar, *forceVar;
        //!read values in a new value and vector
        virtual void readState(STATE s, int rec);
        //!write a new value and vector
        virtual void writeState(STATE s, double time = -1, int rec = -1);


    protected:
        //! Set all of the netcdf dimensions and variables in the file (for creating or writing new files)
        void SetDimVar();
        //! When reading (or writing to an existing files) load in the pre-existing netcdf info
        void GetDimVar();
        //!number of particles in the model
        int N;
        //!size of the vectors
        int dof;
        //! a variable that can be loaded when a state is read
        double val;
        //! a vector for reading doubles in and out
        vector<double> vec;

        //! flags for whether to save and/or load different data structures
        bool velocity, type,force;
        typedef shared_ptr<simpleModel> STATE;
    };
#endif

