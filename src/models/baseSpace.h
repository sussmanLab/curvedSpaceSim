#ifndef baseSpace_H
#define baseSpace_H

#include "std_include.h"
#include <vector>
#include "pointDataType.h"
#include "cgalIncludesAndTypedefs.h"
#include "noiseSource.h"

/*! \file baseSpace.h"
 Degrees of freedom live in a space (e.g.: Euclidean space, a curved manifold, a cube with periodic boundary conditions, etc.)
This base class promises an implementation of
"displaceParticle", which takes a point and and a displacement vector, and updates the position of the point
"distance", which takes a source point, a vector of target points, and fills vectors of distances and tangent vectors at the start and end of the paths.
Because of the eventual target of simulations on curved surfaces, we use the CGAL data types here.
For convenience, baseSpaces also implement the ability to choose a random location (not necessarily uniformly sampled) and a random tangent vector to a point
*/

class baseSpace
    {
    public:
        virtual ~baseSpace() = default;
        //!displace the position of a particle
        virtual void displaceParticle(meshPosition &pos, vector3 &displacementVector) = 0;

        //!move a particle and, if necessary, update the associated velocity vector
        virtual void transportParticleAndVelocity(meshPosition &pos, vector3 &v, vector3 &displacementVector) = 0;

        //!Given a source particle and a vector of target points, determine the geodesic distance and store the start and end path tangents along the paths. The tangents are stored as NORMALIZED vectors
        virtual void distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold) = 0;

        //given a vector of meshPositions, extract some double3 information
        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result)=0;
        //given a vector of meshPositions, extract some double3 information
        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result)=0;

	virtual double getArea()=0; 

        //!Some spaces know that the associated model's default meshPosition structure is already euclidean, but others (e.g., meshes) are not
        bool positionsAreEuclidean = true;

        virtual void randomPosition(meshPosition &p, noiseSource &noise)=0;
        virtual void randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise)=0;
    };
#endif
