#include "functionUtilities.h"
#include "dataTypes.h"

double vectorMagnitude(vector3 v)
    {
    return sqrt(v.squared_length());
    };

vector3 normalize(vector3 v)
    {
    /*
    if (vectorMagnitude(v) == 0) 
        {
	printf("Cannot normalize vector of magnitude 0.\n");
        return v; 
        }
    */
    return v / vectorMagnitude(v);
    }

point3 rotateAboutAxis(point3 p, std::vector<point3> axis, double angle)
    {
    vector3 axisVector(axis[0],axis[1]);
    double axisNorm = axisVector.squared_length();
    axisVector /= sqrt(axisNorm);

    //shift so that we think about rotating the point relative to an origin of the coordinate system at the axis base
    point3 shiftedPoint = p - vector3({0.,0.,0.},axis[0]);
    double3 rotatedTargetCoordinates;
    //We now implement the rodrigues rotation formula, but pre-simplified
    double dotProduct = axisVector[0]*shiftedPoint[0]+axisVector[1]*shiftedPoint[1]+axisVector[2]*shiftedPoint[2];
    double s = sin(angle);
    double c = cos(angle);
    rotatedTargetCoordinates.x = axisVector[0]*dotProduct*(1.-c) + shiftedPoint[0]*c + (axisVector[1]*shiftedPoint[2] - axisVector[2]*shiftedPoint[1])*s;
    rotatedTargetCoordinates.y = axisVector[1]*dotProduct*(1.-c) + shiftedPoint[1]*c + (axisVector[2]*shiftedPoint[0] - axisVector[0]*shiftedPoint[2])*s;
    rotatedTargetCoordinates.z = axisVector[2]*dotProduct*(1.-c) + shiftedPoint[2]*c + (axisVector[0]*shiftedPoint[1] - axisVector[1]*shiftedPoint[0])*s;

    point3 rotatedPoint(rotatedTargetCoordinates.x,rotatedTargetCoordinates.y,rotatedTargetCoordinates.z);
    //shift back to the real coordinate system
    rotatedPoint += vector3({0.,0.,0.},axis[0]);
    return rotatedPoint;
    }

double angleBetween(vector3 v1, vector3 v2)
    {
    v1 = normalize(v1);
    v2 = normalize(v2);
    double oneDotTwo = CGAL::scalar_product(v1,v2);
    return acos(oneDotTwo);
    }

int wrap(int x,int m)
    {
    int ans = x;
    if(x >= m)
        ans = x % m;
    while(ans <0)
        ans += m;
    return ans;
    }

//!fit integers into non-negative domains
int3 wrap(int3 x,int3 m)
    {
    int3 ans;
    ans.x = wrap(x.x,m.x);
    ans.y = wrap(x.y,m.y);
    ans.z = wrap(x.z,m.z);
    return ans;
    }
