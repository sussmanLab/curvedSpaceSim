#ifndef closedMeshSpace_H
#define closedMeshSpace_H

#include "triangulatedMeshSpace.h"

/*!
A closedMeshSpace uses the base triangulatedMeshSpace functions for particle
and vector transport. A closed mesh has no boundaries, so no special functions
are needed
*/
class closedMeshSpace : public triangulatedMeshSpace
    {
    public:
        closedMeshSpace(){};
    };

#endif

