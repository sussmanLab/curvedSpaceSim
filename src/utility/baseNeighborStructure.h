#ifndef baseNeighborStructure_H
#define baseNeighborStructure_H

#include "pointDataType.h"
#include <vector>

/*! \file baseNeighborStructure.h */

/*!
Neighbor structures are those that models can use to populate vectors of vectors of particle indices recording the neighbor j of particle i. This base class provides basic functionality that simple returns all indices other than the target particle index as a potential neighbor (i.e., it does all-to-all neighbors). It should be noted that these classes will typically only return *candidate* neighbors of the particles. For instance, in the case of finding "candidate" neighbors on a curved surface, all particles within a given Euclidean distance might be returned, which a mesh class will later refine. Again, this base class returns all particles as candidate neighbors.
*/
class baseNeighborStructure
    {
    public:
        baseNeighborStructure(){};
        virtual ~baseNeighborStructure() = default;

        virtual void setInteractionRange(double range){maximumInteractionRange = range;};
        //!initialize any helper structure with a set of points (e.g., populating a cell list with a vector of points to sort)
        virtual void initialize(std::vector<meshPosition> &_particles);

        //!populate the candidate lists. Possibly return the largest distance from the source to any target in the candidate list. Include an offset (minIdx) for mpi runs
        virtual double constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset = 0);

        //! A list of indices that will populate the vector of neighbors of particle i
        std::vector<int> indices;
        //! A vector of meshPositions that will populate a list of candidateParticles that *might* be neighbors of a given particle.
        std::vector<meshPosition> particles;

        //!Some neighbor structures, on initialization, need to know whether to require euclidean positions (e.g., cell lists that partition space)
        bool requireEuclideanPositions = false;

        double maximumInteractionRange;
    };
#endif
