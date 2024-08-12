#include "std_include.h"
#include <tclap/CmdLine.h>


#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "noseHooverNVT.h"
#include "velocityVerletNVE.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "simpleModelDatabase.h"
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

using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in curved space!",' ',"V0.9");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<int> saveFrequencyArg("s","saveFrequency","how often a file gets updated",false,100,"int",cmd);
    ValueArg<int> chainLengthArg("l", "M", "length of N-H chain", false, 2, "int", cmd); 
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("e","dt","timestep size",false,.01,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,.2,"double",cmd);
    

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg dangerousSwitch("d","dangerousMeshes","meshes where submeshes are dangerous", cmd, false);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);
    SwitchArg tangentialSwitch("c", "tangentialBCs", "use tangential boundary conditions for open surfaces", cmd, false);

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int maximumIterations = iterationsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    int M = chainLengthArg.getValue(); 
    string meshName = meshSwitchArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    double temperature = temperatureArg.getValue();
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = dangerousSwitch.getValue(); //not used right now
    bool tangentialBCs = tangentialSwitch.getValue(); 

    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
    if(programBranch >0)
        meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,dangerous);
    meshSpace->useTangentialBCs = tangentialBCs; 

    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    configuration->setVerbose(verbose);
    configuration->setSpace(meshSpace);
    cout << "Area: " << meshSpace->getArea() << endl;

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
    if(programBranch >= 1)
        configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly in a small space. Similarly, set random velocities in the tangent plane
    noiseSource noise(reproducible);
    configuration->setRandomParticlePositions(noise);
    configuration->setMaxwellBoltzmannVelocities(noise,temperature);

    //shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,.5);
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);

    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);

    //for now we fix timescale (tau) at 1.0
    cout << "number of N-H masses: " << M << endl; 
    shared_ptr<noseHooverNVT> NVTUpdater=make_shared<noseHooverNVT>(dt, temperature, 1.0, M);
    NVTUpdater->setModel(configuration);
    simulator->addUpdater(NVTUpdater,configuration);
    
    profiler timer("various parts of the code");

    //by default, the simpleModelDatabase will save euclidean positions, mesh positions (barycentric + faceIdx), and particle velocities. See constructor for saving forces and/or particle types as well
    simpleModelDatabase saveState(N,"./testModelDatabase.nc",NcFile::Replace);
    saveState.writeState(configuration,0.0);
    cout << "starting temp: " << NVTUpdater->getTemperatureFromKE() << endl; 
    
    ofstream temperatureFile("nvtTemperatures.csv"); 
 
    double area = meshSpace->getArea(); 
    double density = N/area; 
    double rhoSigma2 = density*maximumInteractionRange*maximumInteractionRange; 
    
    int breakpoint = maximumIterations/100;  
    double sumBPoR = 0;
    int counter = 1; 

    for (int ii = 0; ii < 2*maximumIterations; ++ii)
        {
        timer.start();
        simulator->performTimestep();
        timer.end();
        if(ii%saveFrequency == saveFrequency-1)
            {
            /*
            getFlatVectorOfPositions(configuration,posToSave);
            vvdat.writeState(posToSave,dt*ii);
            */
            saveState.writeState(configuration,dt*ii);
            double nowTemp = NVTUpdater->getTemperatureFromKE();
            printf("step %i T %f \n",ii,nowTemp);
            temperatureFile << ii << ", " << nowTemp << "\n";
	    }
	if ((ii > maximumIterations) && (ii%breakpoint == 0))
       	    {
            double nowTemp = NVTUpdater->getTemperatureFromKE();
	    vector<double> stress(9, 0.0); 
	    simulator->computeMonodisperseStress(stress);
	    //stress trace is pressure, kb = 1, get temperature from average KE, and density is number density above
            sumBPoR += (stress[0]+stress[4]+stress[8])/(density*nowTemp);
	    counter+=1;  
	    }
        }
    
    temperatureFile.close(); 
    timer.print();
    
    double betaPoverRho = sumBPoR/counter; // mean beta*P/rho "after equilibrium" (hopefully)  
    
    cout << "rho sigma squared: " << rhoSigma2 << endl; 
    cout << "beta P over Rho: " << betaPoverRho << endl;
    cout << "{" << rhoSigma2 << ", " << betaPoverRho << "}" << endl;

    string distancesFilename = "../sphere_data/NVTrun_N" + to_string(N) + "_a"+to_string(maximumInteractionRange) + "_distances.csv";
    std::ofstream distancesFile(distancesFilename);

    vector<double> tempDistList;
    vector<vector3> startTangents;
    vector<vector3> endTangents;
    vector<meshPosition> positions = configuration->positions;
    //the below is necessary so that we don't accidentally cut off distances
    double largeDouble = 10000.0;
    meshSpace->setNewSubmeshCutoff(largeDouble);
    meshSpace->useSubmeshingRoutines(false); 
     
    for (int i = 0; i < (positions.size()-1); i++)
        {
	//size -1 because we can't find the ``last'' particle's distance with itself
	vector<meshPosition>::const_iterator first = positions.begin() + i + 1;
        vector<meshPosition>::const_iterator last = positions.end(); 
	vector<meshPosition> targets(first, last);	
	tempDistList.reserve(positions.size() - (i+1)); //largely just readability, a resize call happens within distance 
        meshSpace->distance(positions[i], targets, tempDistList, startTangents, endTangents);
	//as usual, we only care about distances, so we could just discard startTangents and endTangents
        for (double dist: tempDistList)
       	    { 
            distancesFile << dist << ", "; 
	    }
        }
    distancesFile.close(); 

    return 0;
    };
