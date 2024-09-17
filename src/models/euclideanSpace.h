#ifndef euclideanSpace_H
#define euclideanSpace_H

#include "baseSpace.h"

/*! \file euclideanSpace.h" */

/*!
The simplest version of a space, distances are just computed as the difference between points in R^3, and particle displacements are trivial vector addition.
A useful space for debugging, e.g., the force and equation of motion parts of the code base
*/
class euclideanSpace : public baseSpace
    {
    public:
        //!Simple vector addition of particle displacements
        virtual void displaceParticle(meshPosition &pos, vector3 &displacementVector);
        //! in euclideanSpace, velocity vectors are unchanged by particle shifts...this just calls displace particle
        virtual void transportParticleAndVelocity(meshPosition &pos, vector3 &v, vector3 &displacementVector);

        //!Simple vector subtraction for particle separation vectors
        virtual void distance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> &distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent, double distanceThreshold = VERYLARGEDOUBLE);

        //!Ignore the faceIndex field of meshPositions to convert to euclidean positions
        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<double3> &result);
        //!Ignore the faceIndex field of meshPositions to convert from euclidean positions
        virtual void meshPositionToEuclideanLocation(std::vector<meshPosition> &p1, std::vector<meshPosition> &result);
        
        //! just grab a random number in [-.5,.5]^3 for a random position
        virtual void randomPosition(meshPosition &p, noiseSource &noise);
        //!fills v with a random vector (zero mean unit gaussian in each direction)
        virtual void randomVectorAtPosition(meshPosition &p, vector3 &v, noiseSource &noise);
        
        //!returns 0, but could be changed to something more sensible later TODO
        virtual double getArea(); 
    };
#endif
