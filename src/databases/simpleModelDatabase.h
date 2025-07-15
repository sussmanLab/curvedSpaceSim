#ifndef simpleModelDatabase_h
#define simpleModelDatabase_h

#include "baseHDF5Database.h"
#include "simpleModel.h"

/*! \file simpleModelDatabase.h */
//! A database to save states of simpleModels
/*!
There is one unlimited dimension, and each record can store a time and vector information associated with the model
state. The database always saves position information; flags on the constructor control whether other information is
stored
*/
class simpleModelDatabase : public baseHDF5Database
    {
public:
    simpleModelDatabase(int numberOfParticles,
                        string fn = "temp.h5",
                        fileMode::Enum _mode = fileMode::readonly,
                        bool saveVelocities = true,
                        bool saveTypes = false,
                        bool saveForces = false);

    //! Write the current state of the system to the database. If the default value of "rec=-1" is used, just append the
    //! current state to a new record at the end of the database
    virtual void writeState(STATE c, double time = -1.0, int rec = -1);
    //! Read the "rec"th entry of the database into SPV2D state c. If geometry=true, the local geometry of cells
    //! computed (so that further simulations can be run); set to false if you just want to load and analyze
    //! configuration data.
    virtual void readState(STATE c, int rec);

    //! The number of frames that have been saved so far
    unsigned long currentNumberOfRecords();

private:
    typedef shared_ptr<simpleModel> STATE;
    int N; //!< number of points
    int Current; //!< keeps track of the current record when in write mode

    //! a vector of length 1
    std::vector<double> timeVector;
    //! a vector of length 3*N
    std::vector<double> coordinateVector;
    //! a vector of length N
    std::vector<double> doubleVector;
    //! a vector of length N
    std::vector<int> intVector;
    void registerDatasets();

    //! flags for whether to save and/or load different data structures
    bool velocity, type, force;
    };
#endif
