#ifndef geometryCentralMeshSpace_H
#define geometryCentralMeshSpace_H

#include "functionUtilities.h"
#include "meshUtilities.h"
#include "baseSpace.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/surface_centers.h"
#include "geometrycentral/surface/surface_point.h"

/*! \file geometryCentralMeshSpace.h"
* \brief defines an interface to geometry-central mesh-based functionality
*/

/*!
Uses the vector heat method to compute approximate geodesics
*/
class geometryCentralMeshSpace : public baseSpace
    {
    public:
        geometryCentralMeshSpace(){};

        //!load data from an off file and initialize all mesh data structures
        void loadMeshFromFile(std::string filename, bool verbose = false);


        //!Given a particle somewhere on the mesh, displace it in the direction of the vector, wrapping around faces to make it a geodesic displacement
        virtual void displaceParticle(meshPosition &pos, vector3 &displacementVector);
        //!On a mesh, if a particle goes across an edge the velocity vector gets rotated along with it
        virtual void transportParticleAndVelocity(meshPosition &pos, vector3 &v, vector3 &displacementVector);

        //!Given a source particle and a vector of target points, determine the geodesic distance and store the start and end path tangents along the paths
        virtual void distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold= VERYLARGEDOUBLE);

        void surfacePointToEuclideanLocation(geometrycentral::surface::SurfacePoint &p, meshPosition &p1);
        void surfacePointToEuclideanLocation(geometrycentral::surface::SurfacePoint &p, double3 &p1);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result);

        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result);

        virtual void randomPosition(meshPosition &p, noiseSource &noise);

        virtual void randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise);


        void meshPositionToSurfacePoint(meshPosition &p, geometrycentral::surface::SurfacePoint &sp);
        void surfacePointToMeshPosition(geometrycentral::surface::SurfacePoint &sp,meshPosition &p);
        //!given a vector of meshPositions that represent barycentric coordinates, fill a second vector of meshPositions that represent the corresponding R3 positions
        void convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b);

        //data structures
        double3 minVertexPosition;
        double3 maxVertexPosition;
        
        //!GeometryCentral manifold and vertex position data structures
        std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh> mesh;
        std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> geometry;
        std::unique_ptr<geometrycentral::surface::VectorHeatMethodSolver> vectorHeatSolver;
    protected:
        bool verbose = false;
        //data structures associated with submeshing routines
        double maximumDistance = 0;
    };
#endif
