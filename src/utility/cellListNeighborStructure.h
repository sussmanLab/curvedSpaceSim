#ifndef cellListNeighborStructure_H
#define cellListStructure_H

#include  "baseNeighborStructure.h"
#include "hyperRectangularCellList.h"

 /*! \file cellListNeighborStructure.h */

class cellListNeighborStructure : public baseNeighborStructure
    {
    public:
        cellListNeighborStructure(){};
        cellListNeighborStructure(std::vector<double> &minPos, std::vector<double> &maxPos, double gridSize);

        virtual void setInteractionRange(double range);
        virtual void initialize(std::vector<meshPosition> &_particles);

        virtual void constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset = 0);

    protected:
        hyperRectangularCellList cellList;

    };
#endif
