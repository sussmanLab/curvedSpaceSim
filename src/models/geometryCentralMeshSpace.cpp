/*! \file geometryCentralMeshSpace.cpp */
#include "geometryCentralMeshSpace.h"
#include <stdexcept>

void geometryCentralMeshSpace::loadMeshFromFile(std::string filename, bool _verbose)
    {
    verbose = _verbose;
    positionsAreEuclidean = false;
    
    if(verbose)
        {
        printf("loading from file %s\n",filename.c_str());
        }
    std::tie(mesh, geometry) = geometrycentral::surface::readManifoldSurfaceMesh(filename);
    /*
    if(!CGAL::IO::read_polygon_mesh(filename, surface))
        {
        if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, surface))
            {
            std::cerr << "Invalid input file." << std::endl;
            throw std::exception();
            }
        };
    if(!CGAL::is_triangle_mesh(surface))
        {
        std::cerr << "Non-triangular mesh" << std::endl;
        throw std::exception();
        };
    int nFaces = surface.number_of_faces();
    int nVertices = surface.number_of_vertices();
    if(verbose)
        {
        printf("input mesh has %i faces and %i vertices\n",nFaces,nVertices);
        };
    //set domain in which surface lives
    minVertexPosition.x = 0;minVertexPosition.y = 0;minVertexPosition.z = 0;
    maxVertexPosition.x = 0;maxVertexPosition.y = 0;maxVertexPosition.z = 0;
    for (vertexIndex v : surface.vertices())
        {
        point3 p = surface.point(v);
        if(p[0] < minVertexPosition.x)
            minVertexPosition.x = p[0];
        if(p[1] < minVertexPosition.y)
            minVertexPosition.y = p[1];
        if(p[2] < minVertexPosition.z)
            minVertexPosition.z = p[2];
        if(p[0] > maxVertexPosition.x)
            maxVertexPosition.x = p[0];
        if(p[1] > maxVertexPosition.y)
            maxVertexPosition.y = p[1];
        if(p[2] > maxVertexPosition.z)
            maxVertexPosition.z = p[2];
        };
    if(verbose)
        printf("mesh spans (%f,%f,%f) to (%f,%f,%f)\n", minVertexPosition.x,minVertexPosition.y,minVertexPosition.z, maxVertexPosition.x,maxVertexPosition.y,maxVertexPosition.z);

    globalSMSP = make_shared<surfaceMeshShortestPath>(surface);
    globalSMSP->build_aabb_tree(globalTree);
    */
    };

void geometryCentralMeshSpace::meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result)
    {
    }

void geometryCentralMeshSpace::meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result)
    {
    }

void geometryCentralMeshSpace::randomPosition(meshPosition &p, noiseSource &noise)
    {
    }

void geometryCentralMeshSpace::randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise)
    {
    }

void geometryCentralMeshSpace::convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b)
    {
    }


/*
As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex)
*/
void geometryCentralMeshSpace::distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold)
    {
    int nTargets = p2.size();
    distances.resize(nTargets);
    startPathTangent.resize(nTargets);
    endPathTangent.resize(nTargets);

    };

/*
As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex)
*/
void geometryCentralMeshSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
    {
    };

//Code copy-pastes a lot of the displace particle routine above...eventually refactor more nicely
void geometryCentralMeshSpace::transportParticleAndVelocity(meshPosition &pos, vector3 &velocityVector, vector3 &displacementVector)
    {
    };

