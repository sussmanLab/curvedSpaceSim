#ifndef cellListNeighborStructure_H
#define cellListNeighborStructure_H

#include  "baseNeighborStructure.h"
#include "hyperRectangularCellList.h"

 /*! \file cellListNeighborStructure.h */

class cellListNeighborStructure : public baseNeighborStructure
    {
    public:
        cellListNeighborStructure(){};
        //for convenience, construct via vectors or double3s to set the domain size
        cellListNeighborStructure(double3 minPos, double3 maxPos, double gridSize);
        cellListNeighborStructure(std::vector<double> &minPos, std::vector<double> &maxPos, double gridSize);

        virtual void setInteractionRange(double range);
        virtual void initialize(std::vector<meshPosition> &_particles);

        virtual double constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset = 0);

    protected:
        hyperRectangularCellList cellList;

        void findCandidateCellList(meshPosition &p, std::vector<int> &cells);
        double squaredInteractionRange;
    };
#endif
