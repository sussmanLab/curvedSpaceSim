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
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("t","dt","timestep size",false,.01,"double",cmd);
    ValueArg<bool> verboseArg("v", "verbose", "verbosity", false, true, "bool", cmd); 
    ValueArg<bool> submeshArg("s", "submeshed", "whether or not to use submesh", false, true, "bool", cmd); 
    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    
    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int numRemeshings = remeshingsArg.getValue();
    int N = particlesArg.getValue();
    int simIterations = simIterationsArg.getValue();  
    string meshName = meshSwitchArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool submeshed = submeshArg.getValue();
    bool dangerous = false; //not used right now
   
    bool verbose = true; // just always be verbose for tests 

    std::cout << "command line arguments parsed" << std::endl;
    //use above arguments to set up the initial space for testing -- we will later 
    //set whether we should submesh or not
    shared_ptr<triangulatedMeshSpace> meshSpace = make_shared<triangulatedMeshSpace>(); 

    meshSpace->loadMeshFromFile(meshName,verbose);

    //for testing, particles initialized randomly. redo at the start of each check requires restarting the configuration 
    noiseSource noise(reproducible);
    std::cout << "Starting with complexity test ";
    if (submeshed) std::cout << "with submeshing.";
    else std::cout << "without submeshing."; 
    std::cout << std::endl;
    
    //with all-to-all forces    
    meshSpace->useSubmeshingRoutines(submeshed); //indicates whether or not we submesh while generating performance data    
    shared_ptr<simulation> simulator = make_shared<simulation>();
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    simulator->addForce(pairwiseForce);
    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
   
    double startingEdgeLength = meanEdgeLength(meshSpace->surface);
    double targetEdgeLength;
    
    std::ofstream complexity_times("cost_v_complexity.csv");

    for (int i = numRemeshings; i > 0; i--) 
        {
	std::cout << "\nremeshing number: " << numRemeshings + 1 - i << std::endl;
	shared_ptr<simpleModel> configuration = make_shared<simpleModel>(N);
        configuration->setSpace(meshSpace);
        configuration->setRandomParticlePositions(noise);

	double minMaxDim = pow((double)N,(1./3.));
        std::vector<double> minPos(3,-minMaxDim);
        std::vector<double> maxPos(3,minMaxDim);
        if(submeshed)
            {
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
	//remesh at the end so we get the first configuration
        targetEdgeLength = i*startingEdgeLength/numRemeshings; 	
        //meshSpace->isotropicallyRemeshSurface(targetEdgeLength);
	std::string writtenMeshName = "remesh_" + std::to_string(numRemeshings+1-i)+".off";	
        //CGAL::IO::write_OFF(writtenMeshName, meshSpace->surface);

	meshSpace->loadMeshFromFile(writtenMeshName, verbose); 

        }

    };
