#include "cellListNeighborStructure.h"

cellListNeighborStructure::cellListNeighborStructure(std::vector<double> &minPos, std::vector<double> &maxPos, double gridSize)
    {
    cellList.setDomainAndGridSize(minPos,maxPos,gridSize);
    };

void cellListNeighborStructure::setInteractionRange(double range)
    {
    maximumInteractionRange = range;
    squaredInteractionRange = maximumInteractionRange*maximumInteractionRange;
    cellList.setGridSize(range);
    };

void cellListNeighborStructure::initialize(std::vector<meshPosition> &_particles)
    {
    particles = _particles;
    indices.resize(particles.size());
    for(int ii = 0; ii < particles.size(); ++ii)
        {
        indices[ii]=ii;
        }
    cellList.sort(particles);
    };

void cellListNeighborStructure::constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset)
    {
    int neighborNumberGuess = cellList.nMax*2;
    candidateNeighborIndices.clear();
    candidateNeighborIndices.reserve(neighborNumberGuess);
    candidateParticles.clear();
    candidateParticles.reserve(neighborNumberGuess);


    int primaryCellIndex = cellList.positionToCellIndex(p);
    std::vector<int> cellsToSearch;
    cellList.getCellNeighbors(primaryCellIndex, cellsToSearch);


    int cellNumber = cellsToSearch.size();
    for (int cc = 0; cc < cellNumber; ++cc)
        {
        int currentCellIndex = cellsToSearch[cc];
        int localNumberOfCells = cellList.elementsPerCell[currentCellIndex];
        for(int ii = 0; ii < localNumberOfCells; ++ii)
            {
            int idx = cellList.indices[cellList.cellListIndexer(ii,currentCellIndex)];
            if( (particleIndex + offset) != idx)
                {
                if(CGAL::squared_distance(p.x, particles[idx].x) < squaredInteractionRange)
                    {
                    candidateParticles.push_back(particles[idx]);
                    candidateNeighborIndices.push_back(indices[idx]);
                    }
                }
            };
        };
    };
