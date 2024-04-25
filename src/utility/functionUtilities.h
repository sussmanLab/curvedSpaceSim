#ifndef functionUtilities_h
#define functionUtilities_h
#include "dataTypes.h"

#include "cgalIncludesAndTypedefs.h"
//!returns v.v
double vectorMagnitude(vector3 v);

//!rotate a point about an axis (defined by a vector of 2 point3s) by an angle
point3 rotateAboutAxis(point3 p, std::vector<point3> axis, double angle);
vector3 normalize(vector3 v);
double angleBetween(vector3 v1, vector3 v2);

//!fit integers into non-negative domains
int wrap(int x,int m);
//!fit integers into non-negative domains
int3 wrap(int3 x,int3 m);

#endif
