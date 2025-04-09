#include "baseDatabase.h"
#include "std_include.h"
#include <tclap/CmdLine.h>

#include "profiler.h"
#include "noiseSource.h"
#include "euclideanSpace.h"
#include "mpiSimulation.h"
#include "gradientDescent.h"
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
/*
For testing, this version of the code looks only at the simple euclidean
(flat space) parts of the code on a non-mpi model
*/
using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in flat space!",' ',"V0.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<int> saveFrequencyArg("s","saveFrequency","how often a file gets updated",false,100,"int",cmd);
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("t","dt","timestep size",false,.01,"double",cmd);

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int maximumIterations = iterationsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = false; //not used right now

    shared_ptr<euclideanSpace> R3Space=make_shared<euclideanSpace>();
    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    if(verbose)
        configuration->setVerbose(true);
    configuration->setSpace(R3Space);


    //testing cellListNeighborStructure in euclidean, non-periodic spaces...make a cell list with the following minimum and maximum dimensions, and unit grid size
    double minMaxDim = pow((double)N,(1./3.));
    std::vector<double> minPos(3,-minMaxDim);
    std::vector<double> maxPos(3,minMaxDim);
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(minPos,maxPos,maximumInteractionRange);
    configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly in a small space
    noiseSource noise(reproducible);
    configuration->setRandomParticlePositions(noise);

    //shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,.5);
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
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

    valueVectorDatabase vvdat("./flatSpaceTestTrajectory.h5",posToSave.size(),fileMode::replace);
    vvdat.writeState(0,posToSave);

    for (int ii = 0; ii < maximumIterations; ++ii)
        {
        timer.start();
        simulator->performTimestep();
        timer.end();
        if(ii%saveFrequency == saveFrequency-1)
            {
            getFlatVectorOfPositions(configuration,posToSave);
            vvdat.writeState(dt*ii,posToSave);
            double fNorm,fMax;
            fNorm = energyMinimizer->getForceNorm();
            fMax = energyMinimizer->getMaxForce();
            printf("step %i fN %f fM %f\n",ii,fNorm,fMax);
            }
        }

    timer.print();

    return 0;
    };
