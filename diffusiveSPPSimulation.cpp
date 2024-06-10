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
#include "diffusiveSPP.h"

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
    CmdLine cmd("diffusive self-propelled particles!",' ',"V0.1");
    
    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,2,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<int> saveFrequencyArg("s","saveFrequency","how often a file gets updated",false,100,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("e","dt","timestep size",false,.01,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,0.0,"double",cmd);
    ValueArg<double> velocityArg("v", "v0", "self-propulsion velocity", false, 1.0, "double", cmd); 
    ValueArg<double> mobilityArg("q", "mu", "particle mobility", false, 1.0, "double", cmd); 
    ValueArg<double> noiseStrengthArg("u", "noise", "rotational noise strength", false, 0.0, "double", cmd); 
    ValueArg<double> DrArg("D", "Dr", "rotational diffusion rate", false, 1.0, "double", cmd); 
    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg verboseSwitch("w","verbose","output more things to screen ", cmd, false);
    SwitchArg dangerousSwitch("d", "dangerousMeshes", "meshes were submeshes are dangerous", cmd, false);

   

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int maximumIterations = iterationsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    double temperature = temperatureArg.getValue();
    double velocity = velocityArg.getValue();
    double mobility = mobilityArg.getValue();
    double noiseStrength = noiseStrengthArg.getValue();
    double Dr = DrArg.getValue(); 
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = dangerousSwitch.getValue(); //not used right now
    
    int PecletNum = floor(velocity/(Dr*maximumInteractionRange/2)); //peclet number for later naming
    /*
    //simple two test positions for a minimal check/sim -- intended for torus_isotropic_remesh.off surface
    int f1 = 181;
    point3 barys1(0.0338081, 0.873334, 0.0928577);
    meshPosition loc1(barys1,f1);
    int f2 =  945;
    point3 barys2(0.500217, 0.339123, 0.16066);
    meshPosition loc2(barys2,f2);
    vector<meshPosition> testPositions;
    testPositions.reserve(2);
    testPositions.push_back(loc1);
    testPositions.push_back(loc2);
    */

    //set space as a triangulatedMeshSpace
    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
    if(programBranch >0)
        meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,dangerous);

    //set model as "configuration" 
    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    configuration->setVerbose(verbose);
    configuration->setSpace(meshSpace);

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
    if(programBranch >= 1)
        configuration->setNeighborStructure(cellList);

    //below initializes particles randomly, but you can use test positions above if you so choose. 
    //configuration->setParticlePositions(testPositions); //if we defined the two test positions above
    shared_ptr<noiseSource> noise = make_shared<noiseSource>(reproducible); 
    configuration->setRandomParticlePositions(*noise);

    //create force
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);

    //create simulation wrapper
    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);
    
    //define updater -- for diffusive spp, sequence is (double _dt, shared_ptr<noiseSource> noise, double _v0 = 1.0, double _noiseStrength = 1.0, double _mobility = 1.0, double _temperature = 0, bool _reproducible=false)
    shared_ptr<diffusiveSPP> selfPropelledUpdater = make_shared<diffusiveSPP>(dt, noise, velocity, noiseStrength, mobility, temperature, reproducible);
   
    selfPropelledUpdater->setModel(configuration);
    selfPropelledUpdater->initializeRandomVelocities();
    selfPropelledUpdater->setRotationalDiffusionConstant(Dr);  
    
    simulator->addUpdater(selfPropelledUpdater,configuration); 
    
    profiler timer("various parts of the code");

    vector<double> posToSave;
    getFlatVectorOfPositions(configuration,posToSave);
    
    string trajectoryFilename = "./SPP_N"+to_string(N)+"_Pe"+to_string(PecletNum)+".nc";
    vectorValueDatabase vvdat(posToSave.size(),trajectoryFilename,NcFile::Replace);
    vvdat.writeState(posToSave,0);

    for (int ii = 0; ii < maximumIterations; ++ii)
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
                double energyState;
                energyState = pairwiseForce->computeEnergy();
		printf("step %i E %.4g\n ",ii, energyState);
                }
            else
                printf("step %i \n",ii);
            }
	}

    timer.print();

    return 0;
    };
