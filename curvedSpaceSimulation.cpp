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
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "cellListNeighborStructure.h"

void getFlatVectorOfPositions(shared_ptr<simpleModel> model, vector<double> &pos)
    {
    int N = model->N;
    pos.resize(3*N);
    for (int ii = 0; ii < N; ++ii)
        {
        point3 p = model->positions[ii].x;
        pos[3*ii+0] = p[0];
        pos[3*ii+1] = p[1];
        pos[3*ii+2] = p[2];
        }
    };

using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in curved space!",' ',"V0.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int maximumIterations = iterationsArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double dt = 0.01;
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();

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
        point3 p(noise.getRealUniform(-.5,.5),noise.getRealUniform(-.5,.5),noise.getRealUniform(-.5,.5));
        pos[ii].x=p;
        pos[ii].faceIndex=ii;
        if(verbose)
            cout << p[0] <<"  " << p[1] << "  " << p[2] << endl;
        }
    configuration->setParticlePositions(pos);


    //shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,.5);
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,1.0);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);

    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);


    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
    energyMinimizer->setModel(configuration);

    simulator->addUpdater(energyMinimizer,configuration);


    profiler timer("various parts of the code");
vector3 vv;

    vector<double> posToSave;
    getFlatVectorOfPositions(configuration,posToSave);

    vectorValueDatabase vvdat(posToSave.size(),"./testTrajectory.nc",NcFile::Replace);
    vvdat.writeState(posToSave,0);

    for (int ii = 0; ii < maximumIterations/2; ++ii)
        {
        timer.start();
        simulator->performTimestep();
        timer.end();
        if(ii%100 == 99)
            {
            getFlatVectorOfPositions(configuration,posToSave);
            vvdat.writeState(posToSave,dt*ii);
            double fNorm,fMax;
            fNorm = energyMinimizer->squaredTotalForceNorm;
            fMax = energyMinimizer->maximumForceNorm;
            printf("step %i fN %f fM %f\n",ii,fNorm,fMax);
            }
        };

    //testing cellListNeighborStructure in euclidean, non-periodic spaces...make a cell list with the following minimum and maximum dimensions, and unit grid size
    double minMaxDim = pow((double)N,(1./3.));
    std::vector<double> minPos(3,-minMaxDim);
    std::vector<double> maxPos(3,minMaxDim);
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(minPos,maxPos,1.0);
    //note that you can set a new neighbor structure from the beginning, or add one in the middle of the simulation, etc

    configuration->setNeighborStructure(cellList);
    profiler timer2("various parts of the code 2");

    for (int ii = maximumIterations/2; ii < maximumIterations; ++ii)
        {
        timer2.start();
        simulator->performTimestep();
        timer2.end();
        if(ii%100 == 99)
            {
            getFlatVectorOfPositions(configuration,posToSave);
            vvdat.writeState(posToSave,dt*ii);
            double fNorm,fMax;
            fNorm = energyMinimizer->squaredTotalForceNorm;
            fMax = energyMinimizer->maximumForceNorm;
            printf("step %i fN %f fM %f\n",ii,fNorm,fMax);
            }
        };

    timer.print();
    timer2.print();

    return 0;
    };
