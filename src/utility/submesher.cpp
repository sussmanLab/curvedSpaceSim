#include "submesher.h"
#include "std_include.h"

triangleMesh submesher::constructSubmeshFromFaceSet(triangleMesh &mesh, std::set<faceIndex> &faces,std::map<vertexIndex,int> &vertexMap,std::map<faceIndex,int> &faceMap)
    {
    triangleMesh submesh;
    std::set<vertexIndex> vertexIndicesToAdd;
    std::vector<vertexIndex> vidx;
    //get the unique set of vertices
    for(faceIndex currentFace : faces)
        {
        getVertexIndicesFromFace(mesh,currentFace,vidx);
        vertexIndicesToAdd.insert(vidx[0]);
        vertexIndicesToAdd.insert(vidx[1]);
        vertexIndicesToAdd.insert(vidx[2]);
        }

    //estimate the size of the submesh we'll make
    int nVertices = vertexIndicesToAdd.size();
    int nFaces = faces.size();
    int nEdgeEstimate = nFaces*3;
    submesh.reserve(nVertices,nEdgeEstimate,nFaces);

    //get a mapping between 0-indexed vertices we'll add and the corresponding vertex index ordering in the faces of the full mesh. Add the vertices whie we're at it
    int ii = 0;
    for (vertexIndex idx : vertexIndicesToAdd)
        {
        vertexMap.insert(std::make_pair(idx, ii));
        ii+=1;
        point3 vertexPosition = mesh.point(idx);
        submesh.add_vertex(vertexPosition);
        }
    ii = 0;
    for(faceIndex currentFace : faces) 
        {
        getVertexIndicesFromFace(mesh,currentFace,vidx);
        submesh.add_face((vertexIndex) vertexMap[vidx[0]],
                         (vertexIndex) vertexMap[vidx[1]],
                         (vertexIndex) vertexMap[vidx[2]]);
        faceMap.insert(std::make_pair(currentFace,ii));
        ii+=1;
        };

    /* //debugging
    printf("submesh with %i faces and %i vertices; (%i, %i)\n",submesh.number_of_vertices(),submesh.number_of_faces(), nVertices,nFaces);
    std::string outName("./subMesh.off");
    CGAL::IO::write_polygon_mesh(outName, submesh);
    */
    return submesh;
    };

triangleMesh submesher::constructSubmeshFromSourceAndTargets(triangleMesh &mesh, faceLocation &source, std::vector<faceLocation> &targets, double &maximumDistanceFromSource,std::map<vertexIndex,int> &vertexMap, std::map<faceIndex,int> &faceMap)
    {
    faceMap.clear();
    vertexMap.clear();
    double squaredDistanceThreshold = maximumDistanceFromSource*maximumDistanceFromSource;
    //details of the source point
    point3 sourcePoint = PMP::construct_point(source,mesh);
    faceIndex sourceFace = source.first;

    //define some sets to keep track of unique faces that have been visited or should be
    std::set<faceIndex> goalFaces;
    std::set<faceIndex> visitedFaces;
    visitedFaces.insert(sourceFace);

    //first, add to the list of goal faces all faces which contain a target point (and aren't the sourceFace). return early if all of the targets are  in the source face or
    //in the source face and immediately adjacent faces
    for (faceLocation target: targets)
        if (target.first != sourceFace)
            goalFaces.insert(target.first);
    if(goalFaces.empty())
        return constructSubmeshFromFaceSet(mesh,visitedFaces,vertexMap,faceMap);

    std::stack<faceIndex> explorationStack;
    std::vector<point3> faceVertices;
    faceIndex neighboringFace;

    halfedgeIndex hf = mesh.halfedge(sourceFace);
    for (halfedgeIndex hi: halfedges_around_face(hf, mesh))
        {
        neighboringFace = mesh.face(mesh.opposite(hi));
        visitedFaces.insert(neighboringFace);
        explorationStack.push(neighboringFace);
        if(goalFaces.count(neighboringFace)>0)
            goalFaces.erase(neighboringFace);
        };
    if(goalFaces.empty())
        return constructSubmeshFromFaceSet(mesh,visitedFaces,vertexMap,faceMap);

    //Next, move to a breadth-first search of faces that checks every connected face
    faceIndex currentFace;
    while(!explorationStack.empty())
        {
        currentFace = explorationStack.top();
        explorationStack.pop();
        hf = mesh.halfedge(currentFace);
        for(halfedgeIndex hi : halfedges_around_face(hf, mesh))
            {
            neighboringFace = mesh.face(mesh.opposite(hi));
            //avoid forming loops by skipping faces that have already been visited
            if(visitedFaces.count(neighboringFace)>0)
                continue;
            //ignore faces that are completely beyond the cutoff distance
            getVertexPositionsFromFace(mesh,neighboringFace,faceVertices);
            if (CGAL::squared_distance(sourcePoint,faceVertices[0]) > squaredDistanceThreshold &&
                CGAL::squared_distance(sourcePoint,faceVertices[1]) > squaredDistanceThreshold &&
                CGAL::squared_distance(sourcePoint,faceVertices[2]) > squaredDistanceThreshold)
                {
                continue;
                };
            visitedFaces.insert(neighboringFace);
            if(goalFaces.count(neighboringFace)>0)
                {
                goalFaces.erase(neighboringFace);
                }

            explorationStack.push(neighboringFace);
            }
        };
    /*
     The final case is one in which the sphere of maximum distance from the 
     source point intersects part of a triangle but none of its vertices. This 
     can happen for a face at the very boundary, of all visited faces, so we can simply add any remaining goal faces
    */
    for(faceIndex currentFace : goalFaces)
        visitedFaces.insert(currentFace);

    /*//debugging
    if(!goalFaces.empty())
        {
        printf("number of goal faces left: %i\n",goalFaces.size());
        triangleMesh submesh = constructSubmeshFromFaceSet(mesh,visitedFaces,vertexMap,faceMap);
        printf("vertices={");
        for(vertexIndex v : submesh.vertices())
            {
                printPoint(submesh.point(v));printf(",");
            }
        printf("};\n");
        printf("targets={");
        for(int ii = 0; ii < targets.size(); ++ii)
            {
                printPoint(PMP::construct_point(targets[ii],mesh));printf(",");
            }
        printf("};\n");

        printf("goal faces:\n");
        for (faceLocation target: targets)
            printf("%i\t",target.first);
        printf("\n sourceface: %i \n",sourceFace);
        ERRORERROR("exploration stack finished traversal without finding all goal faces. Error");
        }
        */

    return constructSubmeshFromFaceSet(mesh,visitedFaces,vertexMap,faceMap);
    };