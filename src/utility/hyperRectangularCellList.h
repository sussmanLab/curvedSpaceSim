#ifndef hyperRectangularCellList_H
#define hyperRectangularCellList_H

#include <vector>
#include "pointDataType.h"
#include "dataTypes.h"
#include "indexer.h"

/*!
A class that can sort points into a grid of buckets.  This enables local
searches for particle neighbors, etc.  Note that at the moment this class can
only handle hyper-rectangular boxes.  For efficiency, flattens everything and
uses indexers to find the correct place in the various internal data structures
 */
class hyperRectangularCellList
    {
    public:
        hyperRectangularCellList(){};

        //!Given a rectilinear domain and a minimum size of the cells, build up the relevant datastructures
        void setDomainAndGridSize(std::vector<double> &minPos, std::vector<double> &maxPos, double _minimumGridSize);
        //!re-build datastructurs for a given minimum size of the underlying cells
        void setGridSize(double _minimumGridSize);

        //!what are the (up to) 27 cells to search around cellIndex?
        void getCellNeighbors(int cellIndex, std::vector<int> &cellNeighbors);

        //! Loop through the nearby cells of the particles and put an associated index in each cell they belong to
        void sort(std::vector<meshPosition> &p);

        //!What cell index would contain the given position?
        int positionToCellIndex(const meshPosition &p);

        //!A vector containing the number of particles sorted into each cell
        std::vector<unsigned int> elementsPerCell;
        //!A vector containing the particle index in particles of each cell. To be used in conjuction with the indexers below
        std::vector<int> indices;
        //! Indexes the cells in the grid themselves (so the index of the bin corresponding to (i,j,k) is bin = cellIndexer(i,j,k)
        Index3D cellIndexer;
        //!Indexes elements in the cell list (i.e., the third element of cell index i is cellListIndexer(3,i)
        Index2D cellListIndexer;

        //!maximum number of elements in any cell, initialized at some sensible value
        int nMax=6;
        //!do we need to wrap around boundaries when looking for adjacent cells?
        bool periodicSpace = false;
    protected: 


        //! reset the size of various internal data structures
        void resetListSizes();
        //minimum and maximum extents of the domain to be partitioned in each direction
        std::vector<double> minimumPositions;
        std::vector<double> maximumPositions;
        //!total extent in each direction
        double3 domainExtent;
        //! size of the hyperrectangular cells in each direction
        double3 cellSizes;
        //! number of cells in each of the directions
        int3 cellNumbers;

        double currentMinimumGridSize = -1.0;
        int totalCells;
    };
#endif
