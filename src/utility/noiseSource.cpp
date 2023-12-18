#include "noiseSource.h"

/*! \file noiseSource.cpp */

int noiseSource::getInt(int minimum, int maximum)
    {
    int answer;
    uniform_int_distribution<int> uniIntRand(minimum,maximum);
    if (Reproducible)
        answer = uniIntRand(gen);
    else
        answer = uniIntRand(genrd);
    return answer;
    };

double noiseSource::getRealUniform(double minimum, double maximum)
    {
    double answer;
    uniform_real_distribution<double> uniRealRand(minimum,maximum);
    if (Reproducible)
        answer = uniRealRand(gen);
    else
        answer = uniRealRand(genrd);
    return answer;
    };

double noiseSource::getRealNormal(double mean, double std)
    {
    double answer;
    normal_distribution<> normal(mean,std);
    if (Reproducible)
        answer = normal(gen);
    else
        answer = normal(genrd);
    return answer;
    };

void noiseSource::setReproducibleSeed(int _seed)
    {
    RNGSeed = _seed;
    mt19937 Gener(RNGSeed);
    gen = Gener;
#ifdef DEBUGFLAGUP
    mt19937 GenerRd(13377);
    genrd=GenerRd;
#endif
    };
