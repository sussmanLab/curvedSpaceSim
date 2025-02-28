#ifndef diffusiveSPP_H
#define diffusiveSPP_H

#include "baseUpdater.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"

class diffusiveSPP : public updater
    {
    public:
        diffusiveSPP(double _dt, shared_ptr<noiseSource> _noise, double _v0 = 1.0, double _noiseStrength = 1.0, double _mobility = 1.0, double _temperature = 0, bool _reproducible=false)
            {
            deltaT = _dt;
	    v0 = _v0;
	    noiseStrength = _noiseStrength;
	    temperature = _temperature;
	    mobility = _mobility;
            reproducible = _reproducible; 
            noise = _noise;  
	    };
    
	virtual void initializeRandomVelocities(); 
	
	virtual void performUpdate(); 

	virtual void setRotationalDiffusionConstant(double _Dr)
	    {
	    Dr = _Dr;
	    //if uniform: noiseStrength = sqrt(6*Dr*deltaT)/M_PI;
	    noiseStrength = sqrt(2*Dr*deltaT); 
	    cout << "setting noise strength to " << noiseStrength << endl;
	    } 

        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            model->particleShiftsRequireVelocityTransport = true;
            initializeFromModel();
            } 
	
	virtual void setLogFile(shared_ptr<ofstream> file) 
	    {
            logfile = file;
	    }
	
    protected:
        vector<vector3> displacements;	
	shared_ptr<noiseSource> noise;
	shared_ptr<ofstream> logfile; 
	double noiseStrength;
	double v0;
	double mobility;
	double temperature;
	double Dr; 
    };
#endif
