#ifndef functionUtilities_h
#define functionUtilities_h

#include "cgalIncludesAndTypedefs.h"
//!returns v.v
double vectorMagnitude(vector3 v);

//!rotate a point about an axis (defined by a vector of 2 point3s) by an angle
point3 rotateAboutAxis(point3 p, std::vector<point3> axis, double angle);

#endif
