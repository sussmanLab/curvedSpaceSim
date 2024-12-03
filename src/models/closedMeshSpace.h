#ifndef closedMeshSpace_H
#define closedMeshSpace_H

#include "triangulatedMeshSpace.h"

class closedMeshSpace : public closedMeshSpace
    {
    //purely a name guard that clarifies that you're on a surface meant to be closed
    public:
        fullStopMeshSpace(){};
    };

#endif

