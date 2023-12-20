#ifndef pointDataType_H
#define pointDataType_H

/*! \file pointDataType.h
defines a meshPosition class which combines a CGAL Point_3 data type with a face index
*/

#include "cgalIncludesAndTypedefs.h"

class meshPosition
    {
    public:
        meshPosition()
            {
            x = point3(NULL,NULL,NULL);
            faceIndex = -1;
            };
        meshPosition(point3 _x, int _fIdx)
            {
            x = _x;
            faceIndex = _fIdx;
            };


        point3 x;
        int faceIndex;
    protected:
    };
#endif
