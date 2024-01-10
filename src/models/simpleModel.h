#ifndef SIMPLEMODEL_H
#define SIMPLEMODEL_H

#include "std_include.h"
#include "pointDataType.h"
#include "baseSpace.h"
#include "baseNeighborStructure.h"

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

        //!Neighbor structures accelerate the process of finding neighbors...if not set, defaults to all-to-all searches
        virtual void setNeighborStructure(shared_ptr<baseNeighborStructure> _structure)
            {
            neighborStructure = _structure;
            }
        //!move the degrees of freedom
        virtual void moveParticles(vector<vector3> &displacements);
        //!get the number of degrees of freedom, defaulting to the number of cells
        virtual int getNumberOfParticles(){return N;};
        //!some models have internally defined forces. Do everything unusual to compute additional forces... by default, sets forces to zero
        virtual void computeForces(bool zeroOutForces=true);

        virtual void findNeighbors(double maximumInteractionRange);
        virtual void setParticlePositions(vector<meshPosition> &newPositions);

        //!The number of particles
        int N;
        //!particle  positions
        vector<meshPosition> positions;
        //!particle velocities
        vector<vector3> velocities;
        //!Forces on particles
        vector<vector3> forces;

        //!list of list of neighbor indexes
        vector<vector<int>> neighbors;
        //!list of list of separation vectors between neighbors
        vector<vector<vector3>> neighborVectors;
        //!list of list of distances neighbors
        vector<vector<double>> neighborDistances;


        //!particle types
        //GPUArray<int> types;
        //!particle radii
        //GPUArray<scalar> radii;
        //!particle masses
        //GPUArray<scalar> masses;
        
        //sometimes, you just want the actual positions of particles
        vector<double3> euclideanLocations;
        virtual void fillEuclideanLocations();

        void setVerbose(bool v){verbose = v;};
    protected:
        shared_ptr<baseSpace> space;
        shared_ptr<baseNeighborStructure> neighborStructure;
        //!For compatibility across spaces (e.g., mesh-based spaces), a space may have internal functionality to jam euclidean data into a meshPosition. This can be used for neighborStructures
        vector<meshPosition> euclideanMeshPosition;

        bool verbose = false;

    };
typedef shared_ptr<simpleModel> ConfigPtr;
typedef weak_ptr<simpleModel> WeakConfigPtr;
#endif
