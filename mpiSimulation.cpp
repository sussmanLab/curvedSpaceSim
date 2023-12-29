#include "std_include.h"
#include <tclap/CmdLine.h>
#include <mpi.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "euclideanSpace.h"
#include "mpiSimulation.h"
#include "gradientDescent.h"
#include "mpiModel.h"
#include "gaussianRepulsion.h"
#include "vectorValueDatabase.h"

void getFlatVectorOfPositions(shared_ptr<mpiModel> model, vector<double> &pos)
    {
    int N = model->NTotal;
    pos.resize(3*N);
    for (int ii = 0; ii < N; ++ii)
        {
        point3 p = model->globalPositions[ii].x;
        pos[3*ii+0] = p[0];
        pos[3*ii+1] = p[1];
        pos[3*ii+2] = p[2];
        }
    };


/*
In this simplest version of MPI-based parallelization, we note that by far the slowest part of the code is finding neighbors (submeshing, building the source tree).
Thus, we share the entire mesh data structure with all ranks, and we share all particle data with all ranks. Each of R ranks takes charge of roughly 1/R of the particles

 */

using namespace TCLAP;
int main(int argc, char*argv[])
    {
    //First, we handle some MPI initialization:
    int myRank,worldSize;
    int tag=99;
    char message[20];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

    char processorName[MPI_MAX_PROCESSOR_NAME];
    int nameLen;
    MPI_Get_processor_name(processorName, &nameLen);
    int myLocalRank;
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &myLocalRank);

    //then, we set up a basic command line parser with some message and version
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


    if(verbose)
        printf("processes rank %i, local rank %i, world size %i\n",myRank,myLocalRank, worldSize);

    shared_ptr<euclideanSpace> R3Space=make_shared<euclideanSpace>();
    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);

    shared_ptr<mpiModel> configuration=make_shared<mpiModel>(N,myRank,worldSize,verbose);
    configuration->setSpace(R3Space);

    //for testing, just initialize particles randomly in a small space
    noiseSource noise(reproducible);
    vector<meshPosition> pos(N);
    //this block (through "broadcastParticlePositions(pos)") will take the data  on rank 0 and distribute it to all ranks
    for(int ii = 0; ii < N; ++ii)
        {
        point3 p(noise.getRealUniform(),noise.getRealUniform(),noise.getRealUniform());
        pos[ii].x=p;
        pos[ii].faceIndex=ii;
        if(verbose)
            cout << p[0] <<"  " << p[1] << "  " << p[2] << endl;
        }
    configuration->broadcastParticlePositions(pos);


    shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,.5);
    pairwiseForce->setModel(configuration);

    shared_ptr<mpiSimulation> simulator=make_shared<mpiSimulation>(myRank,worldSize);
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);

    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
    energyMinimizer->setModel(configuration);

    simulator->addUpdater(energyMinimizer,configuration);


    profiler timer("various parts of the code");
vector3 vv;

    vector<double> posToSave;
    getFlatVectorOfPositions(configuration,posToSave);
    char dataname[256];
    sprintf(dataname,"parallelTestTrajectory%i.nc",worldSize);
    vectorValueDatabase vvdat(posToSave.size(),dataname,NcFile::Replace);
    vvdat.writeState(posToSave,0);

    for (int ii = 0; ii < maximumIterations; ++ii)
        {
        timer.start();
       // printf("timestep %i on rank %i\n",ii,myRank);
        MPI_Barrier(MPI_COMM_WORLD);
        simulator->performTimestep();
        timer.end();
        if(ii%100 == 99)
            {
            getFlatVectorOfPositions(configuration,posToSave);
            vvdat.writeState(posToSave,dt*ii);
            }
        };
    timer.print();



    MPI_Finalize();

    return 0;
    };
