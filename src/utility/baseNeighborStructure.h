#ifndef baseNeighborStructure_H
#define baseNeighborStructure_H

#include "pointDataType.h"

/*! \file baseNeighborStructure.h */

/*!
Provides basic functionality to find candidate neighbors based on various criteria. In the base version, simply returns all potential indices as candidates
*/

class baseNeighborStructure
    {
    public:
        baseNeighborStructure(){};

        //!populate the candidate lists. Include an offset (minIdx) for mpi runs
        virtual void constructCandidateNeighborList(meshPosition &p, vector<meshPosition> &otherParticles, vector<int> &candidateNeighborIndices, vector<meshPosition> &candidateParticles, int offset = 0);
    };
#endif
