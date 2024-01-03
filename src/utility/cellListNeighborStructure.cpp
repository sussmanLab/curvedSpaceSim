#include "cellListNeighborStructure.h"

cellListNeighborStructure::cellListNeighborStructure(std::vector<double> &minPos, std::vector<double> &maxPos, double gridSize)
    {
    cellList.setDomanAndGridSize(minPos,maxPos,gridSize);
    };

void cellListNeighborStructure::setInteractionRange(double range)
    {
    cellList.setGridSize(gridSize);
    };

void cellListNeighborStructure::initialize(std::vector<meshPosition> &_particles)
    {
    };

void cellListNeighborStructure::constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset = 0)
    {
    };
