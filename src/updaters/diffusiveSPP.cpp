#include "diffusiveSPP.h"
#include "noiseSource.h"
#include "meshUtilities.h"
#include "functionUtilities.h"

//point3 rotateAboutAxis(point3 p, std::vector<point3> axis, double angle)

void diffusiveSPP::initializeRandomVelocities() 
    {
    for (int pp = 0; pp < Ndof; ++pp)
        {
        model->space->randomVectorAtPosition(model->positions[pp], model->velocities[pp], *noise);
        model->velocities[pp] = v0*normalize(model->velocities[pp]);
        }
    };


void diffusiveSPP::performUpdate()
    {
    if(displacements.size() != Ndof)
        displacements.resize(Ndof);
   
    sim->computeForces();
    
    for (int ii = 0; ii < Ndof; ++ii)
        {
	vector3 randomVector(0,0,0); 
	model->space->randomVectorAtPosition(model->positions[ii], randomVector, *noise);   
	displacements[ii] = deltaT*model->velocities[ii] + deltaT*mobility*model->forces[ii] + sqrt(2*deltaT*temperature*mobility)*randomVector; 

	//random angle for rotation from normal distro	
        double rotationAngle = noise->getRealNormal(0, noiseStrength);
	meshPosition iipos = model->positions[ii]; 
        //cout << ii << ", " << iipos.faceIndex << ", pos:" << iipos.x.x() <<", " << iipos.x.y() << ", " << iipos.x.z()  <<  ", vel:" << model->velocities[ii] << endl;

        //past this point, we can offload this to the space's built-in functions  
        vector3 newVelocity = model->velocities[ii];
        vector3 oldVelocity = newVelocity;	
	model->space->rotateVectorAtPosition(model->positions[ii], newVelocity, rotationAngle); 
	model->velocities[ii] = newVelocity;
        
        //cout << ii << ", post rotation velocity " << newVelocity << endl;

	//check to make sure nothing has broken
	double tol = .01; 
	bool breakflag = ((vectorMagnitude(model->velocities[ii]) > (v0 - tol)) && (vectorMagnitude(model->velocities[ii]) < (v0 +tol)));   
	if ( !breakflag ) { 
            //std::vector<point3> result; 
	    //point3 nullPoint(0,0,0); 
	    //result.push_back(nullPoint); 
	    //result.push_back(nullPoint); 
	    //result.push_back(nullPoint);
            //getVertexPositionsFromFace(model->space->surface, iipos.faceIndex, result);
	    cout << "Updater exiting, velocity is " << model->velocities[ii] << "and v0 was " << v0 << endl;
	    cout << "velocity was " << oldVelocity << endl;
	    cout << "face: " << iipos.faceIndex<< endl;
	    //printPoint(result[0]);
	    //printf(", "); 
	    //printPoint(result[1]);
	    //printf(", ");  
	    //printPoint(result[2]);
	    //printf("\n");  
            cout << "particle: " << ii << endl;	    
            cout << "rotation angle: " << rotationAngle << endl;	
	    cout << "velocity magnitude " << vectorMagnitude(model->velocities[ii]) << endl;
	    std::exit(EXIT_FAILURE);
	    } 

        }
    
    sim->moveParticles(displacements);
    };
