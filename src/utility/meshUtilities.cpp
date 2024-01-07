#include "meshUtilities.h"
smspFaceLocation meshPositionToFaceLocation(const meshPosition &p)
    {
    smspBarycentricCoordinates target = {p.x[0],p.x[1],p.x[2]};
    return smspFaceLocation(faceDescriptor(p.faceIndex),target);
    };

void getVertexIndicesFromFace(triangleMesh &mesh, faceIndex i, std::vector<vertexIndex> &result)
    {
    if(result.size()!=3)
        result.resize(3);
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
    if(result.size()!=3)
        result.resize(3);
    halfedgeIndex hf = mesh.halfedge(i);
    int elements=0;
    for(halfedgeIndex hi : halfedges_around_face(hf, mesh))
        {
        vertexIndex vi = source(hi,mesh);
        result[elements] = mesh.point(vi);
        elements +=1;
        };
    };

/*
https://math.stackexchange.com/questions/3995636/point-of-intersection-of-two-lines-in-barycentric-coordinate-system
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

    intersectionScale1 = -(-line2End[1]+line2End[1]*line2Start[0] + line2Start[1] - line2End[0]*line2Start[1])/denominator;
    intersectionScale2 = (-1+line2Start[0]+line2Start[1])/denominator;

    intersectionPoint[0]= line2Start[0] + intersectionScale2*(line2End[0]-line2Start[0]);
    intersectionPoint[1]= line2Start[1] + intersectionScale2*(line2End[1]-line2Start[1]);
    intersectionPoint[2]= line2Start[2] + intersectionScale2*(line2End[2]-line2Start[2]);
    return (intersectionScale1 >= 0 && intersectionScale1 <=1 &&
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

    intersectionPoint[0]= line2Start[0] + intersectionScale2*(line2End[0]-line2Start[0]);
    intersectionPoint[1]= line2Start[1] + intersectionScale2*(line2End[1]-line2Start[1]);
    intersectionPoint[2]= line2Start[2] + intersectionScale2*(line2End[2]-line2Start[2]);
    return (intersectionScale1 >= 0 && intersectionScale1 <=1 &&
            intersectionScale2 >= 0 && intersectionScale2 <=1 );
    };

bool findTriangleEdgeIntersectionInformation(pmpBarycentricCoordinates sourceBarycentricLocation, pmpBarycentricCoordinates targetBarycentricLocation, pmpBarycentricCoordinates &intersectionPoint, std::vector<vertexIndex> vertexList,  std::vector<vertexIndex> &involvedVertex,std::vector<int> &uninvolvedVertex)
    {
    pmpBarycentricCoordinates iCheck;
    bool v1v2Intersection = intersectionBarcentricLinesV1V2(sourceBarycentricLocation,targetBarycentricLocation,iCheck);
    if(v1v2Intersection)
        {
        intersectionPoint = iCheck;
        uninvolvedVertex.push_back(2);
        involvedVertex.push_back(vertexList[0]);
        involvedVertex.push_back(vertexList[1]);
        }
    bool v2v3Intersection = intersectionBarcentricLinesV2V3(sourceBarycentricLocation,targetBarycentricLocation,iCheck);
    if(v2v3Intersection)
        {
        intersectionPoint = iCheck;
        uninvolvedVertex.push_back(0);
        involvedVertex.push_back(vertexList[1]);
        involvedVertex.push_back(vertexList[2]);
        }
    bool v3v1Intersection = intersectionBarcentricLinesV3V1(sourceBarycentricLocation,targetBarycentricLocation,iCheck);
    if(v3v1Intersection)
        {
        intersectionPoint = iCheck;
        uninvolvedVertex.push_back(1);
        involvedVertex.push_back(vertexList[2]);
        involvedVertex.push_back(vertexList[0]);
        }
    }

