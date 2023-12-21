#ifndef SIMPLEMODEL_H
#define SIMPLEMODEL_H

#include "std_include.h"
#include "pointDataType.h"
#include "baseSpace.h"

/*! \file simpleModel.h
 * \brief defines an interface for models that compute forces
 */

//! A base interfacing class that defines common operations
/*!
This provides an interface, guaranteeing that SimpleModel S will provide access to
S.setGPU();
S.getNumberOfParticles();
S.computeForces();
S.moveParticles();
S.returnForces();
S.returnPositions();
S.returnVelocities();
S.returnRadii();
S.returnMasses();
S.returnTypes();
S.spatialSorting();
S.returnAdditionalData();
*/
class simpleModel
    {
    public:
        //!The base constructor requires the number of particles
        simpleModel(int n);
        //!a blank default constructor
        simpleModel(){};
        //!initialize the size of the basic data structure arrays
        void initializeSimpleModel(int n);

        virtual void setSpace(shared_ptr<baseSpace> _space)
            {
            space = _space;
            }
        //!move the degrees of freedom
        virtual void moveParticles(vector<vector3> &displacements);
        //!get the number of degrees of freedom, defaulting to the number of cells
        virtual int getNumberOfParticles(){return N;};
        //!do everything unusual to compute additional forces... by default, sets forces to zero
        virtual void computeForces(bool zeroOutForces=false);

        //void setParticlePositions(vector<dVec> &newPositions);

        //!The number of particles
        int N;
        //!particle  positions
        vector<meshPosition> positions;
        //!particle velocities
        vector<vector3> velocities;
        //!Forces on particles
        vector<vector3> forces;
        //!particle types
        //GPUArray<int> types;
        //!particle radii
        //GPUArray<scalar> radii;
        //!particle masses
        //GPUArray<scalar> masses;

    protected:
        shared_ptr<baseSpace> space;

    };
typedef shared_ptr<simpleModel> ConfigPtr;
typedef weak_ptr<simpleModel> WeakConfigPtr;
#endif
