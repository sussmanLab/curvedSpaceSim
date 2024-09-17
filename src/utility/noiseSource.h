#ifndef noiseSource_H
#define noiseSource_H

#include "std_include.h"
#include "dataTypes.h"

/*! \file noiseSource.h */
//!A class that gives access to a RNG on the cpu
/*!
Provides features to some psuedo-rng functions.  One can call for a random
integer (in a specified range), a random real with a uniform distribution, or a
random real from a normal distribution.
*/
class noiseSource
    {
    public:
        //!base constructor
        noiseSource(bool rep = false)
            {
            Reproducible = rep;
            mt19937 Gener(13377);
        #ifndef DEBUGFLAGUP
            mt19937 GenerRd(rd());
        #else
            mt19937 GenerRd(13377);
        #endif
            gen = Gener;
            genrd=GenerRd;
            }

        //!Get an integer in some range (inclusive of the endpoints)
        int getInt(int minimum, int maximum);
        //!Get a real from uniform distribution in some range.
        double getRealUniform(double minimum =0., double maximum =1.);
        //!Get a real from normal distribution of given mean and standard deviation
        double getRealNormal(double mean =0., double std =1.);

        //!Return a possible barycentric coordinate
        double3 getRandomBarycentricSet();
        //!set reproducibility
        void setReproducible(bool _rep){Reproducible = _rep;};
        //!set the seed on a reproducible RNG run
        void setReproducibleSeed(int _seed);
        //!should the dynamics be reproducible?
        bool Reproducible;
        //!The seed used by the random number generator, when non-reproducible dynamics have been set
        int RNGSeed;
        //!an initializer for non-reproducible random number generation on the cpu
        random_device rd;
        //!A reproducible Mersenne Twister
        mt19937 gen;
        //!A non-reproducible Mersenne Twister
        mt19937 genrd;
    };

#endif
