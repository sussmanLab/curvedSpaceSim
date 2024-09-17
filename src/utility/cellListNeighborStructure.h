#ifndef cellListNeighborStructure_H
#define cellListNeighborStructure_H

#include  "baseNeighborStructure.h"
#include "hyperRectangularCellList.h"

 /*! \file cellListNeighborStructure.h */

/*!
This neighbor structure partitions euclidean 3D space into cells (determined by
the gridSize argument of the constructor) and builds a list of candidate
neighbors by including particles in the 27 cells that include or surround a
given particle's position.  This can also be used as a neighbor structure for
curved space simulations: the euclidean distance used in these cell lists will
be an lower bound on the geodesic distance, so any particles that are within a
given Euclidean distance are a superset of those within a given geodesic
distance.  The meshSpace classes know how to deal with the resulting possibility
of candidate particles on disconnected submeshes.
*/
class cellListNeighborStructure : public baseNeighborStructure
    {
    public:
        cellListNeighborStructure(){};
        //!for convenience, construct via vectors to set the domain size
        cellListNeighborStructure(double3 minPos, double3 maxPos, double gridSize);
        //!..or for convenience, construct via a vector of doubles to set the domain size
        cellListNeighborStructure(std::vector<double> &minPos, std::vector<double> &maxPos, double gridSize);

        virtual void setInteractionRange(double range);
        //!The initialize call prepares data structures and sorts the particles into cells
        virtual void initialize(std::vector<meshPosition> &_particles);

        virtual double constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset = 0);

    protected:
        //!Underlying this class is a cell list that partitions space into cells, and can sort particle positions into these cells
        hyperRectangularCellList cellList;

        //!Given a particle position, return a vector of the 27 cells
        void findCandidateCellList(meshPosition &p, std::vector<int> &cells);
        double squaredInteractionRange;
    };
#endif
