#include "std_include.h"
#include <tclap/CmdLine.h>
#include <mpi.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "mpiSimulation.h"
#include "velocityVerletNVE.h"
#include "mpiModel.h"
#include "harmonicRepulsion.h"
#include "cellListNeighborStructure.h"
#include "vectorValueDatabase.h"

/*
This file was used to generate data of computational time vs particle number, etc, for the paper.

 */

using namespace TCLAP;
int main(int argc, char*argv[])
    {
    //First, we handle some MPI initialization:
    int myRank,worldSize;
    //int tag=99;
    //char message[20];
    //MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

    char processorName[MPI_MAX_PROCESSOR_NAME];
    int nameLen;
    MPI_Get_processor_name(processorName, &nameLen);

    /*
    //local ranks and shmcomm currently not used... perhaps integrate later if warranted
    int myLocalRank;
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &myLocalRank);
    */
    //then, we set up a basic command line parser with some message and version
    CmdLine cmd("Timing runs for curved space simulations!",' ',"V0.9");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,1,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,15,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("e","dt","timestep size",false,.01,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,.2,"double",cmd);

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, false);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int powersOfTwo = particleNumberSwitchArg.getValue();
    int maximumIterations = iterationsArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    double temperature = temperatureArg.getValue();
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();

    vector<double> timingDataSave;
    char dataname[256];
    sprintf(dataname,"timingData_z%i_i%i_a%.2f_mpiRanks%i.nc",programBranch,maximumIterations,maximumInteractionRange,worldSize);


    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
    if(programBranch >=1)
        meshSpace->useSubmeshingRoutines(true,maximumInteractionRange);

    int N;
    for (int i = 4; i < powersOfTwo; ++i)
        {
        N = pow(2,i);
        shared_ptr<mpiModel> configuration=make_shared<mpiModel>(N,myRank,worldSize,verbose);
        configuration->setVerbose(verbose);
        configuration->setSpace(meshSpace);

        //set up the cellListNeighborStructure, which needs to know how large the mesh is
        shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
        if(programBranch >= 1)
            configuration->setNeighborStructure(cellList);

        //just initialize particles randomly on the mesh faces
        noiseSource noise(reproducible);
        configuration->setRandomParticlePositions(noise);
        configuration->setMaxwellBoltzmannVelocities(noise,temperature);

        shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
        pairwiseForce->setModel(configuration);

        shared_ptr<mpiSimulation> simulator=make_shared<mpiSimulation>(myRank,worldSize);
        simulator->setConfiguration(configuration);
        simulator->addForce(pairwiseForce);

        shared_ptr<velocityVerletNVE> nve=make_shared<velocityVerletNVE>(dt);
        nve->setModel(configuration);
        simulator->addUpdater(nve,configuration);

        profiler timer("cost per timestep");

        if(verbose)
            printf("preparing to run for %i steps\n",maximumIterations);
        for (int ii = 0; ii < maximumIterations; ++ii)
            {
            timer.start();
            // printf("timestep %i on rank %i\n",ii,myRank);
            simulator->performTimestep();
            timer.end();
            };

        if(myRank ==0)
            {
            timingDataSave.push_back((double) N);
            timingDataSave.push_back(timer.timing());
            timer.print();
            }

        };
    if(myRank ==0)
        {
        vectorValueDatabase vvdat(timingDataSave.size(),dataname,NcFile::Replace);
        vvdat.writeState(timingDataSave,0);
        };
    MPI_Finalize();

    return 0;
    };
