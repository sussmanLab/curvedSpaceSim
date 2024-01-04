#include "meshUtilities.h"
smspFaceLocation meshPositionToFaceLocation(const meshPosition &p)
    {
    smspBarycentricCoordinates target = {p.x[0],p.x[1],p.x[2]};
    return smspFaceLocation(faceDescriptor(p.faceIndex),target);
    };

void getVertexIndicesFromFace(triangleMesh &mesh, faceIndex i, std::vector<point3> &result)
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
