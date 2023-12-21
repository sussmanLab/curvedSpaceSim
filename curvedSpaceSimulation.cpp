#include "std_include.h"
#include <tclap/CmdLine.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "euclideanSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"

using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in curved space!",' ',"V0.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,50,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double dt = 0.001;
    bool verbose = true;
    bool reproducible = true;

    shared_ptr<euclideanSpace> R3Space=make_shared<euclideanSpace>();
    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);


    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    configuration->setSpace(R3Space);

    //for testing, just initialize particles randomly in a small space
    noiseSource noise(reproducible);
    vector<meshPosition> pos(N);
    for(int ii = 0; ii < N; ++ii)
        {
        point3 p(noise.getRealUniform(),noise.getRealUniform(),noise.getRealUniform());
        pos[ii].x=p;
        pos[ii].faceIndex=0;
        }
    configuration->setParticlePositions(pos);


    shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,1.0);
    pairwiseForce->setModel(configuration);

    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);


    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
    energyMinimizer->setModel(configuration);

    simulator->addUpdater(energyMinimizer,configuration);


    profiler timer("various parts of the code");
vector3 vv;


    timer.start();
    simulator->performTimestep();
    timer.end();
    timer.print();

    return 0;
    };
