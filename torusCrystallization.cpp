#include "std_include.h"
#include <tclap/CmdLine.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "velocityVerletNVE.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "cellListNeighborStructure.h"

#include "pointDataType.h" //for point 3 and related
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

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in curved space!",' ',"V0.9");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> saveFrequencyArg("s","saveFrequency","how often a file gets updated",false,100,"int",cmd); 
    ValueArg<int> toleranceArg("i","crystallizationTolerance","maximum allowed fN exponent",false, 6, "int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<string> positionsFileArg("f", "positionsFile", "filename for xyz positions you want to load", false, "none", "string", cmd); 
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("e","dt","timestep size",false,.01,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,.2,"double",cmd);
    ValueArg<double> energyScaleArg("g", "energyScale", "riesz energy scale", false, 1, "double", cmd);
    ValueArg<double> stiffnessArg("q", "stiffness", "harmonic stiffness", false, 1.0, "double", cmd); 
    
    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);
    SwitchArg dangerousSwitch("d", "dangerousMeshes", "meshes were submeshes are dangerous", cmd, false); 

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int crystallizationTolerance = toleranceArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    string meshName = meshSwitchArg.getValue();
    string positionsFile = positionsFileArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    double temperature = temperatureArg.getValue();
    double energyScale = energyScaleArg.getValue();
    double stiffness = stiffnessArg.getValue();  
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = dangerousSwitch.getValue(); //not used right now

    double maxTolerance = pow(10,-1*crystallizationTolerance);

    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
    if(programBranch >0)
        meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,dangerous);

    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    configuration->setVerbose(verbose);
    configuration->setSpace(meshSpace);

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
    if(programBranch >= 1)
        configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly in a small space. Similarly, set random velocities in the tangent plane
    noiseSource noise(reproducible);
    
    if (positionsFile != "none") {        
       configuration->setMeshPositionsFromR3File(positionsFile, meshSpace->surface);
    }
    else configuration->setRandomParticlePositions(noise);

    configuration->setMaxwellBoltzmannVelocities(noise,temperature);
 
    //riesz potential order is scale, stiffness exponent
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(stiffness,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    //shared_ptr<rieszPotential> pairwiseForce = make_shared<rieszPotential>(energyScale,3.0);
    pairwiseForce->setModel(configuration);

    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);

    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
    shared_ptr<velocityVerletNVE> nve=make_shared<velocityVerletNVE>(dt);
    if(programBranch >=2)
        {
        nve->setModel(configuration);
        simulator->addUpdater(nve,configuration);
        }
    else
        {
        energyMinimizer->setModel(configuration);
        simulator->addUpdater(energyMinimizer,configuration);
        }
    profiler timer("various parts of the code");

    vector<double> posToSave;
    getFlatVectorOfPositions(configuration,posToSave);
    
    string trajectoryFilename = "./torusCrystallization_" + to_string(N) + "_" + to_string(crystallizationTolerance) + ".nc"; 

    vectorValueDatabase vvdat(posToSave.size(),trajectoryFilename,NcFile::Replace);
    vvdat.writeState(posToSave,0);

    double fNorm,fMax, energyState; 
    fNorm = 1;
    int ii = 0;
    while (fNorm > maxTolerance)
        {
        timer.start();
        simulator->performTimestep();
        timer.end();
        if(ii%saveFrequency == saveFrequency-1)
            {
            getFlatVectorOfPositions(configuration,posToSave);
            vvdat.writeState(posToSave,dt*ii);
            if(programBranch <2)
                {
                fNorm = energyMinimizer->getForceNorm();
                fMax = energyMinimizer->getMaxForce();\
		energyState = pairwiseForce->computeEnergy();
                printf("step %i fN %.4g fM %.4g E %.4g\n ",ii,fNorm,fMax, energyState);
                }
            else
                printf("step %i \n",ii);
            }
	ii++; 
        }
    timer.print();    
    std::cout << "stiffness was " << stiffness << std::endl; 
    std::cout << "tolerance exponent was " << crystallizationTolerance << std::endl;
    cout << "final energy was " << pairwiseForce->computeEnergy() << endl;

    string mposFilename = to_string(N) + "_" + to_string(crystallizationTolerance) + "_particle_voronoi_positions.csv";
    std::ofstream meshPosFile(mposFilename); 
    
    configuration->fillEuclideanLocations();
    for (int jj = 0; jj < N; ++jj)
        {
	int fIndex = configuration->positions[jj].faceIndex; 
        double3 p = configuration->euclideanLocations[jj];
	meshPosFile << fIndex << ", " << p.x << ", " << p.y << ", " << p.z << "\n";
        }


    std::ofstream neighborNumFile("neighborNumbers.csv"); 
    configuration->findNeighbors(maximumInteractionRange); 
    vector<int> neighborNumVec;
    neighborNumVec.reserve(configuration->neighbors.size()); 
    for (vector<int> neighbors: configuration->neighbors)
        { 
        neighborNumVec.push_back(neighbors.size()); 
        }  
    std::sort(neighborNumVec.begin(), neighborNumVec.end()); 

    for (int numNeighbors: neighborNumVec) 
        {
        neighborNumFile  <<  numNeighbors << ", ";
        }

    neighborNumFile << "\n" << "Counts:" << "\n"; 
    for (int i = 0; i < 12; i++)
        { 
        neighborNumFile << i << ", " << std::count(neighborNumVec.begin(), neighborNumVec.end(), i) << "\n";	
	}

    return 0;
    };
