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
#include "meshUtilities.h"

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
    std::cout << "parsing command line." << std::endl;
    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("Setting up time vs. number of particles test...",' ',"V0.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> remeshingsArg("i","numRemeshings","number of remeshings to do",false,10,"int",cmd); 
    ValueArg<int> particlesArg("n", "N", "number of particles", false, 100,"int", cmd);
    ValueArg<int> simIterationsArg("j", "simIterations", "number of steps to average over", false, 10, "int", cmd);  
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<string> writeFilePrefix("w", "writeFile", "prefix for written file + meshes", false, "torus", "string", cmd);
    ValueArg<double> deltaTArg("t","dt","timestep size",false,.01,"double",cmd);
    ValueArg<bool> verboseArg("v", "verbose", "verbosity", false, true, "bool", cmd); 
    ValueArg<bool> submeshArg("z", "submeshed", "whether or not to use submesh", false, true, "bool", cmd); 
    ValueArg<bool> newFilesArg("y", "newRemeshFiles", "whether to write new files for remeshes", false, true, "bool", cmd); 
    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
     
    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int numRemeshings = remeshingsArg.getValue();
    int N = particlesArg.getValue();
    int simIterations = simIterationsArg.getValue();  
    string meshName = meshSwitchArg.getValue();
    string writeFile = writeFilePrefix.getValue();
    double dt = deltaTArg.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool submeshed = submeshArg.getValue();
    bool writeNewFiles = newFilesArg.getValue(); 

    bool verbose = true; // just always be verbose for tests 

    std::cout << "command line arguments parsed" << std::endl;
    //use above arguments to set up the initial space for testing -- we will later 
    //set whether we should submesh or not
    shared_ptr<triangulatedMeshSpace> meshSpace = make_shared<triangulatedMeshSpace>(); 

    meshSpace->loadMeshFromFile(meshName,verbose);

    //set max interaction range scaled to the size of the space
    double minMaxDim = pow((double)N,(1./3.));
    std::vector<double> minPos(3,-minMaxDim);
    std::vector<double> maxPos(3,minMaxDim);
    minPos[0] = meshSpace->minVertexPosition.x;
    minPos[1] = meshSpace->minVertexPosition.y;
    minPos[2] = meshSpace->minVertexPosition.z;
    maxPos[0] = meshSpace->maxVertexPosition.x;
    maxPos[1] = meshSpace->maxVertexPosition.y;
    maxPos[2] = meshSpace->maxVertexPosition.z; 
    double diff1 = maxPos[0]-minPos[0];
    double diff2 = maxPos[1]-minPos[1];
    double diff3 = maxPos[2]-minPos[2];
    double spatialExtent = sqrt(diff1*diff1 + diff2*diff2 + diff3*diff3);
    //we scale maximum interaction range roughly with the size of 
    //the mesh such that we get a neighborlist of about the same size every time
    double maximumInteractionRange = spatialExtent/10;
    cout << "max interaction range is " << maximumInteractionRange << endl;

    //for testing, particles initialized randomly. redo at the start of each check 
    //requires restarting the configuration 
    noiseSource noise(reproducible);
    std::cout << "Starting with complexity test ";
    if (submeshed) std::cout << "with submeshing.";
    else std::cout << "without submeshing."; 
    std::cout << std::endl;
    
    //indicate submeshing or not     
    meshSpace->useSubmeshingRoutines(submeshed);      
    //create key simulation primitives 
    shared_ptr<simulation> simulator = make_shared<simulation>();
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);
    simulator->addForce(pairwiseForce);
    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
   
    double startingEdgeLength = meanEdgeLength(meshSpace->surface);
    double targetEdgeLength;
    
    string times_filename = to_string(submeshed) + writeFile + "_cost_v_complexity.csv"; 
    cout << "writing to " << times_filename << endl;
    std::ofstream complexity_times(times_filename);
 
    for (int i = numRemeshings; i > 0; i--) 
        {
	std::cout << "\nremeshing number: " << numRemeshings + 1 - i << std::endl;
	shared_ptr<simpleModel> configuration = make_shared<simpleModel>(N);
        configuration->setSpace(meshSpace);
        configuration->setRandomParticlePositions(noise);
	
	if(submeshed)
            {
	    // had to declare these earlier to set interaction range
            minMaxDim = pow((double)N,(1./3.));
            minPos[0] = meshSpace->minVertexPosition.x;
            minPos[1] = meshSpace->minVertexPosition.y;
            minPos[2] = meshSpace->minVertexPosition.z;
            maxPos[0] = meshSpace->maxVertexPosition.x;
            maxPos[1] = meshSpace->maxVertexPosition.y;
            maxPos[2] = meshSpace->maxVertexPosition.z;
            shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(minPos,maxPos,maximumInteractionRange);
            configuration->setNeighborStructure(cellList);
	    };
        
        pairwiseForce->setModel(configuration);         
	energyMinimizer->setModel(configuration);

        simulator->setConfiguration(configuration);
        simulator->addUpdater(energyMinimizer,configuration);

        profiler timer("remeshed_timer");
        for (int ii = 0; ii < simIterations; ++ii)
                {
		std::cout << "sim step " << ii << std::endl;
                timer.start();
                simulator->performTimestep();
                timer.end();
                //we don't need any trajectory saving because we're just
                //seeing this exact metric
                }
	complexity_times << "\n" << meshSpace->surface.number_of_vertices() << ", " << timer.timing(); 
	//remesh at the end so we get the first configuration -- EDIT HERE TO CHANGE LOWER BOUND
        targetEdgeLength = 3*i*startingEdgeLength/(4*numRemeshings)+startingEdgeLength/4; //minimize at 1/4 of starting edge length to avoid errors from overly complex meshes 	
        if (writeNewFiles) meshSpace->isotropicallyRemeshSurface(targetEdgeLength);
	std::string writtenMeshName = "../remeshes/" + writeFile + "_remesh_" + std::to_string(numRemeshings+1-i)+".off";	
        if (writeNewFiles) CGAL::IO::write_OFF(writtenMeshName, meshSpace->surface);

	meshSpace->loadMeshFromFile(writtenMeshName, verbose); 

        }

    };
