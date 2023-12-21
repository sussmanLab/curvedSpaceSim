#include "euclideanSpace.h"

void euclideanSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
    {
    pos.x += displacementVector;
    };

void euclideanSpace::distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent)
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
        startPathTangent[ii] = Rab;
        endPathTangent[ii] = -Rab;
        };
    };
