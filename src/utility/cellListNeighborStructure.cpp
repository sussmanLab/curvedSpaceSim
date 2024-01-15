#include "cellListNeighborStructure.h"
#include "std_include.h"
cellListNeighborStructure::cellListNeighborStructure(std::vector<double> &minPos, std::vector<double> &maxPos, double gridSize)
    {
    cellList.setDomainAndGridSize(minPos,maxPos,gridSize);
    };

void cellListNeighborStructure::setInteractionRange(double range)
    {
    requireEuclideanPositions = true;
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
//printf("%i %f\n", particles.size(),particles[0].x[0]);
    cellList.sort(particles);
    };

double  cellListNeighborStructure::constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset)
    {
    int neighborNumberGuess = cellList.nMax*2;
    candidateNeighborIndices.clear();
    candidateNeighborIndices.reserve(neighborNumberGuess);
    candidateParticles.clear();
    candidateParticles.reserve(neighborNumberGuess);


    int primaryCellIndex = cellList.positionToCellIndex(p);
    std::vector<int> cellsToSearch;
    cellList.getCellNeighbors(primaryCellIndex, cellsToSearch);

    double maximumDistance = 0;
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
                double dist2 = CGAL::squared_distance(p.x, particles[idx].x);
                if(dist2 < squaredInteractionRange)
                    {
                    candidateParticles.push_back(particles[idx]);
                    candidateNeighborIndices.push_back(indices[idx]);
                    if(dist2 > maximumDistance)
                        maximumDistance = dist2;
                    }
                }
            };
        };
    return maximumDistance;
    };
