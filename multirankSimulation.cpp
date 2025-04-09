#include "std_include.h"
#include <tclap/CmdLine.h>
#include <mpi.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "mpiSimulation.h"
#include "gradientDescent.h"
#include "velocityVerletNVE.h"
#include "mpiModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "cellListNeighborStructure.h"

void getFlatVectorOfPositions(shared_ptr<mpiModel> model, vector<double> &pos)
    {
    int N = model->NTotal;
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
In this simplest version of MPI-based parallelization, we note that by far the slowest part of the code is finding neighbors (submeshing, building the source tree).
Thus, we share the entire mesh data structure with all ranks, and we share all particle data with all ranks. Each of R ranks takes charge of roughly 1/R of the particles

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
    CmdLine cmd("mpi simulations in curved space!",' ',"V0.9");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<int> saveFrequencyArg("s","saveFrequency","how often a file gets updated",false,100,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("e","dt","timestep size",false,.01,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,.2,"double",cmd);

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg dangerousSwitch("d","dangerousMeshes","meshes where submeshes are dangerous", cmd, false);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);

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
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = false; //not used right now


    if(verbose && !dangerous)
        printf("processes rank %i, world size %i\n",myRank, worldSize);
        //printf("processes rank %i, local rank %i, world size %i\n",myRank,myLocalRank, worldSize);

    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
    if(programBranch >=1)
        meshSpace->useSubmeshingRoutines(true,maximumInteractionRange);

    shared_ptr<mpiModel> configuration=make_shared<mpiModel>(N,myRank,worldSize,verbose);
    configuration->setVerbose(verbose);
    configuration->setSpace(meshSpace);

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
    if(programBranch >= 1)
        configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly on the mesh faces
    noiseSource noise(reproducible);
    configuration->setRandomParticlePositions(noise);
    configuration->setMaxwellBoltzmannVelocities(noise,temperature);

    //shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,.5);
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);

    shared_ptr<mpiSimulation> simulator=make_shared<mpiSimulation>(myRank,worldSize);
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
    char dataname[256];
    sprintf(dataname,"parallelTestTrajectory%i.h5",worldSize);
    valueVectorDatabase vvdat(dataname,posToSave.size(),fileMode::replace);
    vvdat.writeState(0,posToSave);
    if(verbose)
        printf("preparing to run for %i steps\n",maximumIterations);
    double fNorm,fMax;
    for (int ii = 0; ii < maximumIterations; ++ii)
        {
        timer.start();
       // printf("timestep %i on rank %i\n",ii,myRank);
        simulator->performTimestep();
        timer.end();
        if(ii%saveFrequency == saveFrequency-1)
            {
            getFlatVectorOfPositions(configuration,posToSave);
            if(myRank ==0)
                vvdat.writeState(dt*ii,posToSave);
            if(programBranch <2)
                {
                fNorm = energyMinimizer->getForceNorm();
                fMax = energyMinimizer->getMaxForce();
                if(myRank ==0)
                    printf("step %i fN %f fM %f\n",ii,fNorm,fMax);
                }
            else
                {
                if(myRank ==0)
                    printf("step %i \n",ii);
                }
            }
        };

    if(myRank ==0)
        {
        timer.print();
        }

    MPI_Finalize();

    return 0;
    };
