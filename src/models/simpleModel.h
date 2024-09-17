#ifndef SIMPLEMODEL_H
#define SIMPLEMODEL_H

#include "std_include.h"
#include "pointDataType.h"
#include "baseSpace.h"
#include "baseNeighborStructure.h"
#include "noiseSource.h"
#include "pointDataType.h" //for point 3 and related
#include "meshUtilities.h" //purely for the r3 position conversions

/*! \file simpleModel.h
 */

//! A base interfacing class that defines common operations
/*!
This provides an interface, guaranteeing that simpleModel S will provide things
like positions, forces, velocities, radii, masses, types, and some additional
data.

Importantly, these models must also implement a moveParticles function in order to interface with the simulation classes, and a findNeighbors function to interface with force classes Other functions exist to help with initialization and allow for operation in both euclidean and mesh spaces
*/
class simpleModel
    {
    public:
        //!The base constructor requires the number of particles
        simpleModel(int n);
        //!a blank default constructor
        simpleModel(){};
        virtual ~simpleModel() = default;
        //!initialize the size of the basic data structure arrays
        void initializeSimpleModel(int n);
        //!set the space pointer to the argument of this function. Can be either a euclidean or mesh space
        virtual void setSpace(shared_ptr<baseSpace> _space)
            {
            space = _space;
            }

        //!Neighbor structures accelerate the process of finding neighbors...if not set, defaults to all-to-all searches
        virtual void setNeighborStructure(shared_ptr<baseNeighborStructure> _structure)
            {
            neighborStructure = _structure;
            }
        //!move the degrees of freedom. Calls the space's function to do this
        virtual void moveParticles(vector<vector3> &displacements);
        //!get the number of degrees of freedom, defaulting to the number of cells
        virtual int getNumberOfParticles(){return N;};
        //!some models have internally defined forces. Do everything unusual to compute additional forces... by default, sets forces to zero
        virtual void computeForces(bool zeroOutForces=true);
        //!Provide a routine to use the neighbor structures to find neighbors within a given distance of partcles
        virtual void findNeighbors(double maximumInteractionRange);
        //! set the particle positions to the argument of this function
        virtual void setParticlePositions(vector<meshPosition> &newPositions);
        //! Avoid numerical precision issues by clamping near-zero barycentric coordinate to zero
    	virtual void clampBarycentricCoordinatesToFace(point3 &barycentricWeights);
        //!takes a vector of euclidean positions, and appends to "simPositions" a set of (barycentric coordinates, face index) of the closest point on the mesh.
        virtual void R3PositionsToMeshPositions(triangleMesh &mesh, vector<point3> r3positions, vector<meshPosition> &simPositions);
        //!A function of convenience to load particle positions from a file
        virtual void setMeshPositionsFromR3File(string filename, triangleMesh &mesh);
        //!uses the space to randomly set particle positions...hence, requires that a space is already set
        virtual void setRandomParticlePositions(noiseSource &noise);
        //!Set particle velocities to be drawn from a maxwell-boltzmann distribution
        virtual void setMaxwellBoltzmannVelocities(noiseSource &noise, double T);

        //!The space the model lives in
        shared_ptr<baseSpace> space;

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
        vector<int> types;
        //!particle radii TODO
        //vector<double> radii;
        //!particle masses
        vector<double> masses;
        
        //sometimes, you just want the actual positions of particles
        vector<double3> euclideanLocations;
        //!fill the euclideanLocations vector with the euclidean locations of particles 
        virtual void fillEuclideanLocations();

        //!optionally turn on some print statements
        void setVerbose(bool v){verbose = v;};
        //!some updaters require parallel transport of velocity vectors
        bool particleShiftsRequireVelocityTransport = false;

        //!tolerance used for barycentricCoordinate clamping    :with
        double clampTolerance = 0.00000000000001;//10^-14 as a current threshold for numerical tolerance. 
    protected:
        //! a pointer to a data structure that can be used to accelerate finding nearby neighbors
        shared_ptr<baseNeighborStructure> neighborStructure;
        //!For compatibility across spaces (e.g., mesh-based spaces), a space may have internal functionality to jam euclidean data into a meshPosition. This can be used for neighborStructures
        vector<meshPosition> euclideanMeshPosition;

        //!optionally turn on some print statements
        bool verbose = false;
    };
typedef shared_ptr<simpleModel> ConfigPtr;
typedef weak_ptr<simpleModel> WeakConfigPtr;
#endif
