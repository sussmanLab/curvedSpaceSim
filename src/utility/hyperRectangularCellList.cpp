#include "hyperRectangularCellList.h"

#include "functions.h"

/*!
This routine currently picks an even integer of cells in each dimension, close to but larger than the desired size, that fit in the box.
 */
void hyperRectangularCellList::setGridSize(double _minimumGridSize)
    {
    if(currentMinimumGridSize == _minimumGridSize)
        return;
    currentMinimumGridSize = _minimumGridSize;

    domainExtent.x = maximumPositions[0] - minimumPositions[0];
    domainExtent.y = maximumPositions[1] - minimumPositions[1];
    domainExtent.z = maximumPositions[2] - minimumPositions[2];

    cellNumbers.x = floor(domainExtent.x / currentMinimumGridSize);
    if(cellNumbers.x%2==1) cellNumbers.x -=1;
    cellSizes.x = domainExtent.x / cellNumbers.x;

    cellNumbers.y = floor(domainExtent.y / currentMinimumGridSize);
    if(cellNumbers.y%2==1) cellNumbers.y -=1;
    cellSizes.y = domainExtent.y / cellNumbers.y;

    cellNumbers.z = floor(domainExtent.z / currentMinimumGridSize);
    if(cellNumbers.z%2==1) cellNumbers.z -=1;
    cellSizes.z = domainExtent.z / cellNumbers.z;

    totalCells = cellNumbers.x*cellNumbers.y*cellNumbers.z;

    cellListIndexer = Index2D(nMax,totalCells);
    cellIndexer = Index3D(cellNumbers.x,cellNumbers.y,cellNumbers.z);

    resetListSizes();
    };

void hyperRectangularCellList::resetListSizes()
    {
    //resize elementsPerCell and set to zero
    elementsPerCell.resize(totalCells);
    for (int ii = 0; ii < totalCells; ++ii)
        elementsPerCell[ii] = 0;
    //resize indices and set to zero
    cellListIndexer = Index2D(nMax,totalCells);
    indices.resize(cellListIndexer.size());
    for(int ii = 0; ii < cellListIndexer.size(); ++ii)
        indices[ii] = 0;
    }

void hyperRectangularCellList::setDomainAndGridSize(std::vector<double> &minPos, std::vector<double> &maxPos, double _minimumGridSize)
    {
    minimumPositions = minPos;
    maximumPositions = maxPos;

    setGridSize(_minimumGridSize);
    };

int hyperRectangularCellList::positionToCellIndex(const meshPosition &p)
    {
    int3 cIdx;
    //use maxes and mins to make sure the cell indexer is in bounds
    cIdx.x = std::max(0, std::min(cellNumbers.x-1, (int) floor((p.x[0]+minimumPositions[0])/cellSizes.x )));
    cIdx.y = std::max(0, std::min(cellNumbers.y-1, (int) floor((p.x[1]+minimumPositions[1])/cellSizes.y )));
    cIdx.z = std::max(0, std::min(cellNumbers.z-1, (int) floor((p.x[2]+minimumPositions[2])/cellSizes.z )));

    return cellIndexer(cIdx);
    };

void hyperRectangularCellList::sort(std::vector<meshPosition> &p)
    {
    //will loop through particles and put them in cells...
    //if there are more than Nmax particles in any cell, will need to recompute.
    int3 bin;
    int binIndex, offset, cellListPos;
    int N = p.size();
    int nmax = nMax;
    bool recompute = true;
    while(recompute)
        {
        recompute = false;
        //reset elementsPerCell, indices, and cellListIndexer
        resetListSizes();
        for (int ii = 0; ii < N; ++ii)
            {
            binIndex = positionToCellIndex(p[ii]);
            offset = elementsPerCell[binIndex];
            if(offset < nMax && !recompute)
                {
                cellListPos = cellListIndexer(offset,binIndex);
                indices[cellListPos] = ii;
                }
            else
                {
                nmax = std::max(nMax,offset+1);
                nMax = nmax;
                recompute = true;
                }
            elementsPerCell[binIndex] += 1;
            };
        }
    cellListIndexer = Index2D(nMax,totalCells);
    };

//return a vector of the cells to search if you're interested in cellIndex
void hyperRectangularCellList::getCellNeighbors(int cellIndex, std::vector<int> &cellNeighbors)
    {
    cellNeighbors.reserve(27);
    int3 location = cellIndexer.inverseIndex(cellIndex);
    int xMin, xMax, yMin, yMax, zMin, zMax;
    if(periodicSpace)
        {
        xMin = std::max(0,location.x-1);
        xMax = std::min(cellNumbers.x-1,location.x+1);
        yMin = std::max(0,location.y-1);
        yMax = std::min(cellNumbers.y-1,location.y+1);
        zMin = std::max(0,location.z-1);
        zMax = std::min(cellNumbers.z-1,location.z+1);
        for(int xx = xMin; xx < xMax; ++xx)
            for(int yy = yMin; yy < yMax; ++yy)
                for(int zz = zMin; zz < zMax; ++zz)
                    cellNeighbors.push_back(cellIndexer(xx,yy,zz));
        }
    else
        {
        xMin = location.x-1;
        xMax = location.x+1;
        yMin = location.y-1;
        yMax = location.y+1;
        zMin = location.z-1;
        zMax = location.z+1;
        for(int xx = xMin; xx < xMax; ++xx)
            for(int yy = yMin; yy < yMax; ++yy)
                for(int zz = zMin; zz < zMax; ++zz)
                    cellNeighbors.push_back(cellIndexer(
                                                    wrap(xx,cellNumbers.x),
                                                    wrap(xx,cellNumbers.y),
                                                    wrap(xx,cellNumbers.z)));
        };
    };
