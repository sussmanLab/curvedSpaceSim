#include "std_include.h"
#include <tclap/CmdLine.h>
#include <iostream>
#include <fstream>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "cellListNeighborStructure.h"

void getFlatVectorOfPositions(shared_ptr<simpleModel> model, vector<double> &pos)
    {
    int N = model->N;
    model->fillEuclideanLocations();
    pos.resize(3*N);
    for (int ii = 0; ii < N; ++ii)
        {
        double3 p = model->euclideanLocations[ii];
        pos[3*ii+0] = p.x;
        pos[3*ii+1] = p.y;
        pos[3*ii+2] = p.z;
        }
    };

using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("Setting up time vs. number of particles test...",' ',"V1.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> iterationsArg("i","sim_iterations","number of performTimestep calls to make",false,5,"int",cmd); // default to 5 as we'll do this over and over
    ValueArg<int> doublingsArg("n", "particle_doublings", "number of times to double particle number", false, 10, "int", cmd); 
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<double> interactionRangeArg("a","interactionRange","range of the interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("t","dt","timestep size",false,.01,"double",cmd);

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int maximumIterations = iterationsArg.getValue();
    int doublings = doublingsArg.getValue(); 
    string meshName = meshSwitchArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    bool reproducible = reproducibleSwitch.getValue();
   
    bool verbose = true; // just always be verbose for tests 

    std::cout << "command line arguments parsed" << std::endl;
    //use above arguments to set up the initial space for testing -- we will later 
    //set whether we should submesh or not
    shared_ptr<triangulatedMeshSpace> meshSpace = make_shared<triangulatedMeshSpace>(); 

    meshSpace->loadMeshFromFile(meshName,verbose);

    //for testing, particles initialized randomly. redo at the start of each check requires restarting the configuration 
    noiseSource noise(reproducible);
    std::cout << "Starting tests" << std::endl;
   
    shared_ptr<simulation> simulator = make_shared<simulation>();
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    simulator->addForce(pairwiseForce);
    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
 
    //multiplying by 100 guarantees we don't miss numbers after the decimal in the file name  
    int printRange = 100*maximumInteractionRange;
    string filename = "cost_v_n" + to_string(printRange) + ".csv";
    ofstream tvn_times(filename);

     
    //comment /* here to prevent all-to-all test

    //with all-to-all forces    
    meshSpace->useSubmeshingRoutines(false);    
    /*
    shared_ptr<simulation> simulator = make_shared<simulation>();
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    simulator->addForce(pairwiseForce);
    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
    */  
    
    tvn_times << "All-to-All"; 
    for (int i = 1; i <= doublings; i++) 
        {
	int N = pow(2,i); 
        //at the moment, i have put everything that depends on N here -- 
	//assuming we only use the cellList structure when submeshing is 
	//used
	shared_ptr<simpleModel> configuration = make_shared<simpleModel>(N);
	configuration->setSpace(meshSpace);
	configuration->setRandomParticlePositions(noise);
	
        pairwiseForce->setModel(configuration);
    	energyMinimizer->setModel(configuration);

	simulator->setConfiguration(configuration);
	simulator->clearUpdaters(); 
    	simulator->addUpdater(energyMinimizer,configuration);        

        profiler aa_timer("all to all"); 
	for (int ii = 0; ii < maximumIterations; ++ii)
        	{	
        	aa_timer.start();
        	simulator->performTimestep();
        	aa_timer.end();
		//we don't need any trajectory saving because we're just 
		//seeing this exact metric 
        	}
        tvn_times << "\n" << N << ", " << aa_timer.timing(); 	
	aa_timer.print(); // this will spit out the average time per timestep  
        	
	}
    //stop comment here to prevent all to all test

    tvn_times << "\nSubmeshed, Fixed Interaction Range";

    //with submeshed forces     
    meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,false);

    //preliminaries for celllist so we don't recompute them every step 
        
    for (int i = 1; i <= doublings; i++) 
        {
        int N = pow(2,i); 
        //because we change particle number every doubling, everything that depends on N goes here 
        shared_ptr<simpleModel> sm_configuration = make_shared<simpleModel>(N);
	sm_configuration->setSpace(meshSpace);
	sm_configuration->setRandomParticlePositions(noise);
		
	double minMaxDim = pow((double)N,(1./3.));
        std::vector<double> minPos(3,-minMaxDim);
        std::vector<double> maxPos(3,minMaxDim);
        minPos[0] = meshSpace->minVertexPosition.x;
        minPos[1] = meshSpace->minVertexPosition.y;
        minPos[2] = meshSpace->minVertexPosition.z;
        maxPos[0] = meshSpace->maxVertexPosition.x;
        maxPos[1] = meshSpace->maxVertexPosition.y;
        maxPos[2] = meshSpace->maxVertexPosition.z;
        shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(minPos,maxPos,maximumInteractionRange);
	
	sm_configuration->setRandomParticlePositions(noise);
        sm_configuration->setNeighborStructure(cellList);     

                pairwiseForce->setModel(sm_configuration);
    	energyMinimizer->setModel(sm_configuration);

	simulator->setConfiguration(sm_configuration);
        simulator->clearUpdaters(); 
    	simulator->addUpdater(energyMinimizer,sm_configuration);


        profiler sm_timer("submeshed"); 
	for (int ii = 0; ii < maximumIterations; ++ii)
        	{	
        	sm_timer.start();
        	simulator->performTimestep();
        	sm_timer.end();
		//we don't need any trajectory saving because we're just 
		//seeing this exact metric 
        	}
        tvn_times << "\n" << N << ", " << sm_timer.timing(); 
	sm_timer.print(); // this will spit out the average time per timestep  	
	} 

    return 0;
    };
