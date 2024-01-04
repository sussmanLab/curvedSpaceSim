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
