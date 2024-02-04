/*! \file geometryCentralMeshSpace.cpp */
#include "geometryCentralMeshSpace.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include <stdexcept>

using namespace geometrycentral;
using namespace geometrycentral::surface;

void geometryCentralMeshSpace::loadMeshFromFile(std::string filename, bool _verbose)
    {
    verbose = _verbose;
    positionsAreEuclidean = false;
    
    if(verbose)
        {
        printf("loading from file %s\n",filename.c_str());
        }
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filename);

    geometry->requireVertexPositions();
    geometry->requireVertexIndices();
    geometry->requireFaceNormals();
    geometry->requireFaceIndices();
    geometry->requireFaceTangentBasis();
    geometry->requireHalfedgeVectorsInVertex();
    geometry->requireHalfedgeVectorsInFace();

    vectorHeatSolver = std::make_unique<VectorHeatMethodSolver>(*geometry,1.0); 

    int nFaces = mesh->nFaces();
    int nVertices = mesh->nVertices();
    if(verbose)
        {
        printf("input mesh has %i faces and %i vertices\n",nFaces,nVertices);
        };
    //set domain in which surface lives
    minVertexPosition.x = 0;minVertexPosition.y = 0;minVertexPosition.z = 0;
    maxVertexPosition.x = 0;maxVertexPosition.y = 0;maxVertexPosition.z = 0;
    //for (Vertex v : geometry->vertices())
    for (int ii = 0; ii < nVertices; ++ii)
        {
        Vector3 p = geometry->vertexPositions[ii];
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
        printf("mesh spans (%f,%f,%f) to (%f,%f,%f)\n", minVertexPosition.x,
                        minVertexPosition.y,minVertexPosition.z, maxVertexPosition.x,
                        maxVertexPosition.y,maxVertexPosition.z);
    };

void geometryCentralMeshSpace::surfacePointToEuclideanLocation(SurfacePoint &p, meshPosition &p1)
    {
    Face f = p.face;
    Vector3 baryCoords = p.faceCoords;
    Vector3 ans{0.,0.,0.};
    int ii = 0;
    Halfedge he = f.halfedge();
    Vector3 pA = geometry->vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = geometry->vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = geometry->vertexPositions[he.vertex()];
    ans = baryCoords[0]*pA+baryCoords[1]*pB+baryCoords[2]*pC;
    p1.x = point3(ans[0], ans[1], ans[2]);
    p1.faceIndex = f.getIndex();
    };

void geometryCentralMeshSpace::surfacePointToEuclideanLocation(SurfacePoint &p, double3 &p1)
    {
    Face f = p.face;
    Vector3 baryCoords = p.faceCoords;
    Vector3 ans{0.,0.,0.};
    int ii = 0;
    Halfedge he = f.halfedge();
    Vector3 pA = geometry->vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = geometry->vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = geometry->vertexPositions[he.vertex()];
    ans = baryCoords[0]*pA+baryCoords[1]*pB+baryCoords[2]*pC;
    p1.x = ans[0];
    p1.y = ans[1];
    p1.z = ans[2];
    };

void geometryCentralMeshSpace::meshPositionToSurfacePoint(meshPosition &p, geometrycentral::surface::SurfacePoint &sp)
    {
    Vector3 baryCoords{p.x[0],p.x[1],p.x[2] };
    Face f = mesh->face(p.faceIndex);
    sp = SurfacePoint(f,baryCoords);
    };

void geometryCentralMeshSpace::surfacePointToMeshPosition(geometrycentral::surface::SurfacePoint &sp,meshPosition &p)
    {
    p.x = point3(sp.faceCoords[0],sp.faceCoords[1],sp.faceCoords[2]);
    p.faceIndex = sp.face.getIndex();
    };

void geometryCentralMeshSpace::meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result)
    {
    if(result.size()!=p1.size())
        result.resize(p1.size());
    SurfacePoint sp;
    for(int ii = 0; ii < p1.size();++ii)
        {
        meshPositionToSurfacePoint(p1[ii],sp);
        surfacePointToEuclideanLocation(sp,result[ii]);
        };
    };

void geometryCentralMeshSpace::meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result)
    {
    if(result.size()!=p1.size())
        result.resize(p1.size());
    SurfacePoint sp;
    for(int ii = 0; ii < p1.size();++ii)
        {
        meshPositionToSurfacePoint(p1[ii],sp);
        surfacePointToEuclideanLocation(sp,result[ii]);
        };
    };

void geometryCentralMeshSpace::convertToEuclideanPositions(std::vector<meshPosition> &a, std::vector<meshPosition> &b)
    {
    meshPositionToEuclideanLocation(a,b);
    };

void geometryCentralMeshSpace::randomPosition(meshPosition &p, noiseSource &noise)
    {
    double3 baryPoint;
    baryPoint.x=noise.getRealUniform();
    baryPoint.y=noise.getRealUniform(0,1-baryPoint.x);
    baryPoint.z=1-baryPoint.x-baryPoint.y;
    p.x = point3(baryPoint.x,baryPoint.y,baryPoint.z);
    p.faceIndex = noise.getInt(0,mesh->nFaces()-1);
    }

void geometryCentralMeshSpace::randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise)
    {
    //compute the normal to the face containing point p
    SurfacePoint sourceLocation;
    meshPositionToSurfacePoint(p,sourceLocation);
    Vector3 sourceNormal = geometry->faceNormals[sourceLocation.face];
    //construct an orthogonal vector by hand
    Vector3 orthogonalizer;
    if(sourceNormal[0] ==0 && sourceNormal[1] ==0 )
        orthogonalizer = Vector3{sourceNormal[1]-sourceNormal[2], sourceNormal[2]-sourceNormal[0], sourceNormal[0]-sourceNormal[1]};
    else
        orthogonalizer = Vector3{sourceNormal[1],-sourceNormal[0],0};
    //make this orthogonal vector a unit vector
    orthogonalizer = orthogonalizer.normalize();
    //grab the second tangent vector (which should already be a unit vector)
    Vector3 tangent2 = cross(sourceNormal,orthogonalizer);

    Vector3 result = noise.getRealNormal()*orthogonalizer+noise.getRealNormal()*tangent2;

    //finally, choose a gaussian weight of each of the two orthogonal vectors in the face's tangent
    v = vector3(result[0],result[1],result[2]);
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

    SurfacePoint sourcePoint;
    meshPositionToSurfacePoint(p1,sourcePoint);
    VertexData<Vector2> logMap = vectorHeatSolver->computeLogMap(sourcePoint);

double3 test;
surfacePointToEuclideanLocation(sourcePoint, test);
printf("source = {%f,%f,%f};\n",test.x,test.y,test.z);

    Vector2 logMapV;
    for (int ii =0; ii < nTargets; ++ii)
        {
        Vector2 logMapResult{0.,0.};
        SurfacePoint targetPoint;
        meshPositionToSurfacePoint(p2[ii],targetPoint);
surfacePointToEuclideanLocation(targetPoint, test);
printf("target = {%f,%f,%f};\n",test.x,test.y,test.z);
        //identify the vertices associated with the target face
        int internalCoordinate = 0;
        double dist = 0;
        for (Halfedge he : targetPoint.face.adjacentHalfedges())
            {
            //compute value of logmap
            logMapV = logMap[he.vertex()];
            //compute the change of basis to bring it back to the face (?)
            Vector2 rotation = (geometry->halfedgeVectorsInFace[he] / geometry->halfedgeVectorsInVertex[he]).normalize();
printf("rotationVector = {%f,%f};\n",rotation[0],rotation[1]);
            //accumulate the result
            logMapResult += targetPoint.faceCoords[internalCoordinate] * rotation * logMapV;
            //logMapResult += targetPoint.faceCoords[internalCoordinate] *  logMapV;
printf("currentLogMap = {%f,%f};\npartialResult  = {%f,%f};\n",
                    logMapV[0],logMapV[1],(rotation * logMapV)[0],(rotation * logMapV)[1]);
            dist+= norm((targetPoint.faceCoords[internalCoordinate] * logMapV));
            internalCoordinate +=1;
            };
printf("weightedNorm=%f;\n",dist);
        distances[ii] = norm(logMapResult);
        Vector3 basisX = geometry->faceTangentBasis[sourcePoint.face][0];
        Vector3 basisY = geometry->faceTangentBasis[sourcePoint.face][1];
        logMapResult = logMapResult/distances[ii];
        Vector3 tangentAns = logMapResult[0]*basisX+logMapResult[1]*basisY;
        startPathTangent[ii] = vector3(tangentAns[0],tangentAns[1],tangentAns[2]);
        //empty data!!!!!
        endPathTangent[ii] = vector3(0,0,0);
            /*
        Halfedge he = targetPoint.face.halfedge();
        logMapV1 = logMap[he.vertex()];
        he = he.next();
        logMapV2 = logMap[he.vertex()];
        he = he.next();
        logMapV3 = logMap[he.vertex()];
            */

printf("distance = %f;\ntangent={%f,%f,%f};\n",distances[ii], tangentAns[0],tangentAns[1],tangentAns[2]);

        }

    //https://geometry-central.net/surface/algorithms/vector_heat_method/
    };

/*
As a reminder: in the routine below, we assume that the point3 member of all of the meshPositions
(i.e., p1.x), is actually just carrying around the real numbers corresponding to the barycentric
coordinates of that point in the corresponding faceIndex (i.e., p1.faceIndex)
*/
void geometryCentralMeshSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
    {
    SurfacePoint sp;
    meshPositionToSurfacePoint(pos,sp);
    Vector3 displacement{displacementVector[0],displacementVector[1],displacementVector[2]};
    Vector3 basisX = geometry->faceTangentBasis[sp.face][0];
    Vector3 basisY = geometry->faceTangentBasis[sp.face][1];
    Vector2 tangentVectorDisplacement{dot(displacement,basisX),dot(displacement,basisY)};
    SurfacePoint pathEndpoint = traceGeodesic(*geometry, sp, tangentVectorDisplacement).endPoint;

    surfacePointToMeshPosition(pathEndpoint,pos);
    };

//Code copy-pastes a lot of the displace particle routine above...eventually refactor more nicely
void geometryCentralMeshSpace::transportParticleAndVelocity(meshPosition &pos, vector3 &velocityVector, vector3 &displacementVector)
    {
    UNWRITTENCODE("position and velocity transport not implemented yet for geometryCentralMesh.");
    };

