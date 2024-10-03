#include "fireMinimization.h"

void fireMinimization::minimizeByFire()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
    sim->computeForces();
    if (useFWithExclusions)
        {
	vector<int> exclusions; 
        model->findBoundaryParticles(exclusions);
        forceMax = getMaxForceWithExclusions(exclusions);
        }	
    else forceMax = getMaxForce();
    iterations= 0;
    while((iterations < maximumIterations) && forceMax > forceCutoff)
        {
        iterations +=1;
        velocityVerletFirstHalfStep();
        velocityVerletSecondHalfStep();

        fireStep();
	if (useFWithExclusions)
            {
            vector<int> exclusions2;
            model->findBoundaryParticles(exclusions2);
            forceMax = getMaxForceWithExclusions(exclusions2);
            }
        else forceMax = getMaxForce();
        }
    //cout << "maximum force is: " << forceMax << endl;
    };

double fireMinimization::updaterVectorDotProduct(vector<vector3> &v1, vector<vector3> &v2)
    {
    vector<double> dotProd(1);
    dotProd[0] = 0.0;
    for (int ii = 0; ii < v1.size(); ++ii)
        {
        dotProd[0] += v1[ii]*v2[ii];
        }
    //define a lambda which is just addition
    sim->manipulateUpdaterData(dotProd,
                         [](double x, double y)-> double
                                        {
                                        return x+y;
                                        });
    return dotProd[0];
    };

void fireMinimization::fireStep()
    {
    power = 0.0;
    if (useFWithExclusions)
        {
        vector<int> exclusions;
        model->findBoundaryParticles(exclusions);
        forceMax = getForceNormWithExclusions(exclusions);
        }
    else forceNorm = updaterVectorDotProduct(model->forces,model->forces);
    velocityNorm = updaterVectorDotProduct(model->velocities,model->velocities);
    power = updaterVectorDotProduct(model->forces,model->velocities);

    double scaling = 0.0;
    if(forceNorm>0)
        scaling = sqrt(velocityNorm/forceNorm);

    //adjust velocities according to the fire algorithm 
    for (int ii = 0; ii < Ndof; ++ii)
        model->velocities[ii] = (1-alpha)*model->velocities[ii] + alpha*scaling * model->forces[ii];

    //test if we need to reset velocities
    if(power > 0)
        {
        if(nSinceNegativePower > nMin)
            {
            deltaT = min(deltaT*deltaTInc,deltaTMax);
            alpha = alpha*alphaDec;
            alpha = max(alpha,alphaMin);
            };
        nSinceNegativePower +=1;
        }
    else
        {
        nSinceNegativePower = 0;
        deltaT = deltaT*deltaTDec;
        deltaT = max(deltaT,deltaTMin);
        alpha = alphaStart;
        model->velocities.clear();
        model->velocities.resize(Ndof,vector3(0.,0.,0.));
        };
    };

void fireMinimization::setFIREParameters(int _maximumIterations,double _deltaT, double _alphaStart, double _deltaTMax, double _deltaTMin, double _deltaTInc, double _deltaTDec, double _alphaDec, int _nMin, double _forceCutoff, double _alphaMin)
    {
    maximumIterations = _maximumIterations;
    alphaStart = _alphaStart;
    deltaTMax = _deltaTMax;
    deltaTInc = _deltaTInc;
    deltaTMin = _deltaTMin;
    deltaTDec = _deltaTDec;
    alphaDec = _alphaDec;
    forceCutoff = _forceCutoff;
    alphaMin = _alphaMin;
    nMin = _nMin;
    alpha = alphaStart;
    };
