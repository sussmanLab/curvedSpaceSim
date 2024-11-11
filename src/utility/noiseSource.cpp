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
/*!
To be within a triangle, we choose u from [0,1], v from [0,1-u], and then set w = 1-u-v
*/
double3 noiseSource::getRandomBarycentricSet()
    {
    double3 ans;
    double u,v;
    u=getRealUniform(0,1);
    v=getRealUniform(0,1-u);
    ans.x=u;
    ans.y=v;
    ans.z=1-u-v;
    return ans;
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
