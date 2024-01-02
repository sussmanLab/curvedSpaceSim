#ifndef INDEXER_H
#define INDEXER_H

#include "dataTypes.h"

/*! \file indexer.h */

//!Switch between a 2D array to a flattened, 1D index 
class Index2D
    {
    public:
        Index2D(unsigned int w=0) : width(w), height(w) {}
        Index2D(unsigned int w, unsigned int h) : width(w), height(h) {}

        unsigned int operator()(unsigned int i, unsigned int j) const
            {
            return j*width + i;
            }
        //!Return the number of elements that the indexer can index
        unsigned int size() const
            {
            return width*height;
            }

        //!Get the width
        unsigned int getW() const
            {
            return width;
            }

        //!get the height
        unsigned int getH() const
            {
            return height;
            }

        unsigned int width;
        unsigned int height;
    };

//!Switch between a 3D array to a flattened, 1D index 
class Index3D
    {
    public:
        Index3D(unsigned int w=0){setSizes(w, w ,w);};
        Index3D(int w, int h, int d){setSizes(w,h,d);};

        void setSizes(int w, int h, int d)
            {
            sizes.x = w; sizes.y = h; sizes.z = d;
            numberOfElements = sizes.x * sizes.y * sizes.z;
            intermediateSizes.x=1;
            intermediateSizes.y = intermediateSizes.x*sizes.x;
            intermediateSizes.z = intermediateSizes.y*sizes.y;
            };

        unsigned int operator()(const int x, const int y, const int z) const
            {
            return x*intermediateSizes.x + y*intermediateSizes.y + z*intermediateSizes.z;
            };

        unsigned int operator()(const int3 &i) const
            {
            return i.x*intermediateSizes.x + i.y*intermediateSizes.y + i.z*intermediateSizes.z;
            };

        //! the index i corresponds to what int3 triple under the operator Index3D(int3 )?
        int3 inverseIndex(int i)
            {
            int3 ans;
            int z0 = i;
            ans.x = z0%sizes.x;
            z0= (z0-ans.x)/sizes.x;
            ans.y = z0%sizes.y;
            z0=(z0-ans.y)/sizes.y;
            ans.z = z0%sizes.z;
            return ans;
            };

        //!Return the number of elements that the indexer can index
        unsigned int size() const
            {
            return numberOfElements;
            };

        //!Get the iVec of sizes
        int3 getSizes() const
            {
            return sizes;
            };

        int3 sizes;
        int3 intermediateSizes;
        unsigned int numberOfElements;
    };
#endif
