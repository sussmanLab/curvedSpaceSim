#include "meshUtilities.h"
#include "std_include.h"

smspFaceLocation meshPositionToFaceLocation(const meshPosition &p)
    {
    smspBarycentricCoordinates target = {p.x[0],p.x[1],p.x[2]};
    return smspFaceLocation(faceDescriptor(p.faceIndex),target);
    };
/*!
Note that both getVertex... functions assume the results vector is already the correct size)
 */
void getVertexIndicesFromFace(triangleMesh &mesh, const faceIndex &i, std::vector<vertexIndex> &result)
    {
    halfedgeIndex hf = mesh.halfedge(i);
    int elements=0;
    for(halfedgeIndex hi : halfedges_around_face(hf, mesh))
        {
        vertexIndex vi = source(hi,mesh);
        result[elements] = vi;
        elements +=1;
        };
    };
void getVertexPositionsFromFace(triangleMesh &mesh, faceIndex i, std::vector<point3> &result)
    {
    halfedgeIndex hf = mesh.halfedge(i);
    int elements=0;
    for(halfedgeIndex hi : halfedges_around_face(hf, mesh))
        {
        vertexIndex vi = source(hi,mesh);
        result[elements] = mesh.point(vi);
        elements +=1;
        };
    };

double meanEdgeLength(triangleMesh mesh, bool verbose /*=false*/)
    {
    //finds the mean edge length of all the edges in a mesh via an edge iterator

    //initial typedefs for return/printed quantities
    double mean = 0.0, min, max, length;
    double count = 0.0; bool init = true;

    //ranges have iterator objects assocated with their begin and end, but don't themselves have edges
    triangleMesh::Halfedge_range es = mesh.halfedges();
    for (auto eIter = es.begin(); eIter != es.end(); ++eIter)
        {
        //iterator is a pointer, get an index by checking what it points to
        halfedgeIndex e = *eIter;

        point3 a = mesh.point(mesh.source(e)); //source yields the source vertex of an edge as a vertex_descriptor
        point3 b = mesh.point(mesh.target(e));

        length = CGAL::sqrt(CGAL::squared_distance(a, b));
        ++count;
        if (init){
            mean = min = max = length;
            init = false;
        }
        else{
            if (length < min) min = length;
            if (length > max) max = length;
        }
        mean += length;
        }

   mean /= count;

   if (verbose) std::cout << "edge length min, max, mean: " << min << " " << max << " " << mean << "\n";

   return mean;
   };

double triangleArea(point3 v1, point3 v2, point3 v3)
    {
    vector3 side1 = vector3(v1, v2);
    vector3 side2 = vector3(v1, v3);
    return vectorMagnitude(CGAL::cross_product(side1,side2))/2.0;
    };

double meanTriangleArea(triangleMesh mesh)
    {
    double mean = 0.0;
    int count = 0;
    auto facelist = mesh.faces();
    for (faceIndex face: facelist)
        {
        ++count;
        std::vector<point3> vs;
        vs.reserve(3);
        getVertexPositionsFromFace(mesh,face,vs);
        mean += triangleArea(vs[0],vs[1],vs[2]);
        };
    return mean/count;
    };

double totalArea(triangleMesh mesh)
    {
    double area = 0.0;
    auto facelist = mesh.faces();
    for (faceIndex face: facelist)
        {
        std::vector<point3> vs;
        vs.reserve(3);
        getVertexPositionsFromFace(mesh,face,vs);
        area += triangleArea(vs[0],vs[1],vs[2]);
        };
    return area;
    };

/*simple clamp function to take barycentric coordinate near the boundary of a face
 *and place it firmly within the face in question. Be careful this is not used when 
 *the barycentric coordinates are/might be very negative, as it can obscure meaningful errors. 
 */
void belowZeroClamp(pmpBarycentricCoordinates &baryPoint, double tol)
    {
    double clampedBarySum = 0;
    for (int i = 0; i < 3; i++)
        {
        baryPoint[i] = max(baryPoint[i], tol);
	clampedBarySum += baryPoint[i]; 
        }
    for (int i = 0; i < 3; i++)
        {
        baryPoint[i] = baryPoint[i]/clampedBarySum;
        }
    }

void nearZeroClamp(pmpBarycentricCoordinates &baryPoint, double tol) 
    {
    double clampedBarySum = 0;
    for (int i = 0; i < 3; i++)
        {
        if (baryPoint[i] > -tol && baryPoint[i] < tol) 
	    {
	    baryPoint[i] = tol;
	    }
	clampedBarySum += baryPoint[i]; 
        }
    for (int i = 0; i < 3; i++)
        {
        baryPoint[i] = baryPoint[i]/clampedBarySum;
        }
    } 

/*
Let a line in barycentric coordinates be 
l(t,start,end) = start + t(end-start),
where start = (u,v,w) for w=1-u-v,
then given two lines, solve for an intersection point between the given segments (i.e., 
0 <= t1 <= 1 and 
0 <= t2 <= 1
*/
bool intersectionOfLinesInBarycentricCoordinates(pmpBarycentricCoordinates line1Start, pmpBarycentricCoordinates line1End, pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint)
    {
    double intersectionScale1, intersectionScale2,denominator;
    denominator = (line1Start[0]-line1End[0])*(line2End[1]-line2Start[1]) - (line2End[0]-line2Start[0])*(line1Start[1]-line1End[1]);
    if(denominator ==0)
        return false;

    intersectionScale1 = (line1Start[0]*(line2End[1]-line2Start[1])+line2Start[0]*(line1Start[1]-line2End[1])+line2End[0]*(line2Start[1]-line1Start[1]))/denominator;
    intersectionScale2 = (line1Start[0]*(line1End[1]-line2Start[1])+line1End[0]*(line2Start[1]-line1Start[1])+line2Start[0]*(line1Start[1]-line1End[1]))/denominator;

    intersectionPoint[0]= line2Start[0] + intersectionScale2*(line2End[0]-line2Start[0]);
    intersectionPoint[1]= line2Start[1] + intersectionScale2*(line2End[1]-line2Start[1]);
    intersectionPoint[2]= line2Start[2] + intersectionScale2*(line2End[2]-line2Start[2]);
    return (intersectionScale1 >= 0 && intersectionScale1 <=1 &&
            intersectionScale2 >= 0 && intersectionScale2 <=1 );
    };

bool intersectionBarycentricLinesV1V2(pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint)
    {
    double intersectionScale1, intersectionScale2,denominator;
    denominator = line2Start[0]+line2Start[1] - line2End[0]-line2End[1];
    if(denominator ==0)
        return false;

    intersectionScale1 = (-line2End[1]+line2End[1]*line2Start[0] + line2Start[1] - line2End[0]*line2Start[1])/denominator;
    intersectionScale2 = (-1+line2Start[0]+line2Start[1])/denominator;
//printf("\n 12 scales: %f\t%f\n",intersectionScale1,intersectionScale2);
    intersectionPoint[0]= line2Start[0] + intersectionScale2*(line2End[0]-line2Start[0]);
    intersectionPoint[1]= line2Start[1] + intersectionScale2*(line2End[1]-line2Start[1]);
    intersectionPoint[2]= line2Start[2] + intersectionScale2*(line2End[2]-line2Start[2]);
    return (intersectionScale1 >=0  && intersectionScale1 <=1 &&
            intersectionScale2 >= 0 && intersectionScale2 <=1 );
    };
bool intersectionBarycentricLinesV2V3(pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint)
    {
    double intersectionScale1, intersectionScale2,denominator;
    denominator = -line2End[0] + line2Start[0];
    if(denominator ==0)
        return false;

    intersectionScale1 = -(line2End[0] - line2Start[0] + line2End[1]*line2Start[0] - line2End[0]*line2Start[1])/denominator;
    intersectionScale2 = (line2Start[0])/denominator;
//printf("\n 23 scales: %f\t%f\n",intersectionScale1,intersectionScale2);

    intersectionPoint[0]= line2Start[0] + intersectionScale2*(line2End[0]-line2Start[0]);
    intersectionPoint[1]= line2Start[1] + intersectionScale2*(line2End[1]-line2Start[1]);
    intersectionPoint[2]= line2Start[2] + intersectionScale2*(line2End[2]-line2Start[2]);
    return (intersectionScale1 >= 0 && intersectionScale1 <=1 &&
            intersectionScale2 >= 0 && intersectionScale2 <=1 );
    };

bool intersectionBarycentricLinesV3V1(pmpBarycentricCoordinates line2Start, pmpBarycentricCoordinates line2End, pmpBarycentricCoordinates &intersectionPoint)
    {
    double intersectionScale1, intersectionScale2,denominator;
    denominator = line2Start[1]-line2End[1];
    if(denominator ==0)
        return false;

    intersectionScale1 = (line2End[0]*line2Start[1] - line2End[1]*line2Start[0])/denominator;
    intersectionScale2 = line2Start[1]/denominator;
//printf("\n 31 scales: %f\t%f\n",intersectionScale1,intersectionScale2);

    intersectionPoint[0]= line2Start[0] + intersectionScale2*(line2End[0]-line2Start[0]);
    intersectionPoint[1]= line2Start[1] + intersectionScale2*(line2End[1]-line2Start[1]);
    intersectionPoint[2]= line2Start[2] + intersectionScale2*(line2End[2]-line2Start[2]);
    return (intersectionScale1 >= 0 && intersectionScale1 <=1 &&
            intersectionScale2 >= 0 && intersectionScale2 <=1 );
    };

bool findTriangleEdgeIntersectionInformation(pmpBarycentricCoordinates sourceBarycentricLocation, pmpBarycentricCoordinates targetBarycentricLocation, pmpBarycentricCoordinates &intersectionPoint, std::vector<vertexIndex> vertexList,  halfedgeIndex previousHalfEdge, triangleMesh &surface, std::vector<vertexIndex> &involvedVertex,std::vector<int> &uninvolvedVertex)
    {
    pmpBarycentricCoordinates iCheck;
    bool answer = false;
    if(!(previousHalfEdge == surface.halfedge(vertexList[0],vertexList[1]) || previousHalfEdge == surface.halfedge(vertexList[1],vertexList[0])))
        {
        bool v1v2Intersection = intersectionBarycentricLinesV1V2(sourceBarycentricLocation,targetBarycentricLocation,iCheck);
        if(v1v2Intersection )
            {
            intersectionPoint = iCheck;
            uninvolvedVertex.push_back(2);
            involvedVertex.push_back(vertexList[0]);
            involvedVertex.push_back(vertexList[1]);
            answer=true;
            }
        };
    if(!(previousHalfEdge == surface.halfedge(vertexList[1],vertexList[2]) || previousHalfEdge == surface.halfedge(vertexList[2],vertexList[1])))
        {
        bool v2v3Intersection = intersectionBarycentricLinesV2V3(sourceBarycentricLocation,targetBarycentricLocation,iCheck);
        if(v2v3Intersection )
            {
            intersectionPoint = iCheck;
            uninvolvedVertex.push_back(0);
            involvedVertex.push_back(vertexList[1]);
            involvedVertex.push_back(vertexList[2]);
            answer=true;
            }
        };
    if(!(previousHalfEdge == surface.halfedge(vertexList[2],vertexList[0]) || previousHalfEdge == surface.halfedge(vertexList[0],vertexList[2])))
        {
        bool v3v1Intersection = intersectionBarycentricLinesV3V1(sourceBarycentricLocation,targetBarycentricLocation,iCheck);
        if(v3v1Intersection )
            {
            intersectionPoint = iCheck;
            uninvolvedVertex.push_back(1);
            involvedVertex.push_back(vertexList[2]);
            involvedVertex.push_back(vertexList[0]);
            answer=true;
            }
        };
    return answer;
    }

/*!
Submeshing is very helpful in these calculations, but one then needs to map
between barycentric coordinates of a particle on a face of the submesh and the
barycentric coordinates of that same particle on the original mesh.  This
function does just that.
*/
void convertBarycentricCoordinates(triangleMesh &mesh1, triangleMesh &mesh2, std::unordered_map<faceIndex,int> &faceMap, smspFaceLocation &locationToConvert)
    {
    std::vector<point3> vertexPositions1(3);
    std::vector<point3> vertexPositions2(3);
    getVertexPositionsFromFace(mesh1,locationToConvert.first,vertexPositions1);
    getVertexPositionsFromFace(mesh2,(faceIndex) faceMap[locationToConvert.first],vertexPositions2);
    std::unordered_map<point3,int> sourcePointVertexOrderMap;
    sourcePointVertexOrderMap.insert(std::make_pair(vertexPositions1[0],0));
    sourcePointVertexOrderMap.insert(std::make_pair(vertexPositions1[1],1));
    sourcePointVertexOrderMap.insert(std::make_pair(vertexPositions1[2],2));
    smspBarycentricCoordinates originalCoordinates = locationToConvert.second;
    smspBarycentricCoordinates newCoordinates;
    newCoordinates[0] =originalCoordinates[sourcePointVertexOrderMap[vertexPositions2[0]]];
    newCoordinates[1] =originalCoordinates[sourcePointVertexOrderMap[vertexPositions2[1]]];
    newCoordinates[2] =originalCoordinates[sourcePointVertexOrderMap[vertexPositions2[2]]];
    locationToConvert.second = newCoordinates;
    locationToConvert.first = (faceIndex)faceMap[locationToConvert.first];
    }

void computePathDistanceAndTangents(surfaceMeshShortestPath *smsp, smspFaceLocation &targetPoint, double &distance, vector3 &startPathTangent, vector3 &endPathTangent)
    {
    std::vector<point3> pathPoints;
    shortestPathResult geodesic = smsp->shortest_path_points_to_source_points(targetPoint.first, targetPoint.second,  std::back_inserter(pathPoints));
    distance = std::get<0>(geodesic);
    //Note that the path goes from the target to source, so if we want to know path tangent at the source for force calculation, we must use the *end* of points[]
    if(distance < 0)
        return;
    int pathSize = pathPoints.size();
    
    //debug statement only
    //cout << "path points:\n" <<endl;
    //for (point3 p: pathPoints) cout << p << "\n" << endl;

    startPathTangent = -vector3(pathPoints[pathSize-2],pathPoints[pathSize-1]);
    endPathTangent = -vector3(pathPoints[0],pathPoints[1]);
    //normalize path tangents
    double normalization = sqrt(startPathTangent.squared_length());
    startPathTangent /= normalization;
    normalization = sqrt(endPathTangent.squared_length());
    endPathTangent /= normalization;
    }


void printPoint(point3 a, bool precise)
    {
    if (precise) printf("{%.16g,%.16g,%.16g}",a[0],a[1],a[2]);
    else printf("{%f,%f,%f}",a[0],a[1],a[2]);
    };

void printBary(smspBarycentricCoordinates a, bool precise)
    {
    if (precise) printf("{%.32g,%.32g,%.32g}",a[0],a[1],a[2]);
    else printf("{%f,%f,%f}",a[0],a[1],a[2]);
    };
