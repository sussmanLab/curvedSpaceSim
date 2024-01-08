#include "euclideanSpace.h"

void euclideanSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
    {
    pos.x += displacementVector;
    };

void euclideanSpace::meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result)
    {
    if(result.size()!=p1.size())
        result.resize(p1.size());
    for(int ii = 0; ii < p1.size();++ii)
        {
        result[ii].x=p1[ii].x[0];
        result[ii].y=p1[ii].x[1];
        result[ii].z=p1[ii].x[2];
        };
    };

void euclideanSpace::distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent)
    {
    int nTargets = p2.size();
    distances.resize(nTargets);
    startPathTangent.resize(nTargets);
    endPathTangent.resize(nTargets);
    for (int ii = 0; ii < nTargets; ++ii)
        {
        vector3 Rab(p1.x, p2[ii].x);//CGAL notation for the vector from p1 to p2[ii]
        double dist = sqrt(Rab.squared_length());
        distances[ii] = dist;
        startPathTangent[ii] = (1./dist)*Rab;
        endPathTangent[ii] = (-1./dist)*Rab;
        };
    };
