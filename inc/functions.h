#ifndef functions_H
#define functions_H

#include "dataTypes.h"

//!fit integers into non-negative domains
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
#endif
