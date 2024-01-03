#ifndef baseNeighborStructure_H
#define baseNeighborStructure_H

#include "pointDataType.h"
#include <vector>

/*! \file baseNeighborStructure.h */

/*!
Provides basic functionality to find candidate neighbors based on various criteria. In the base version, simply returns all indices other than the target particleIndex as potential neighbors (i.e., does all-to-all neighbors
*/

class baseNeighborStructure
    {
    public:
        baseNeighborStructure(){};

        virtual void setInteractionRange(double range){};
        //!initialize any helper structure with a set of points (e.g., populating a cell list with a vector of points to sort)
        virtual void initialize(std::vector<meshPosition> &_particles);

        //!populate the candidate lists. Include an offset (minIdx) for mpi runs
        virtual void constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset = 0);

        std::vector<int> indices;
        std::vector<meshPosition> particles;
    };
#endif
