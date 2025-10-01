#include "std_include.h"
#include <tclap/CmdLine.h>


#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "closedMeshSpace.h"
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
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,1,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<int> saveFrequencyArg("f","saveFrequency","how often a file gets updated",false,100,"int",cmd);
    ValueArg<int> chainLengthArg("l", "M", "length of N-H chain", false, 2, "int", cmd); 
    ValueArg<int> sourcesArg("k", "nSources", "number of particles to collect ditsances from", false, 1, "int", cmd);
    ValueArg<int> sampleNoArg("c","sampleNo","sample identifier for this set of parameters",false,1,"int",cmd); 
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/sphere_radius1.off","string",cmd);
    ValueArg<string> savePathArg("s","fileSavePath","path to where data will be saved",false,"./","string",cmd);
    ValueArg<double> areaFractionArg("a","areaFraction","percent of mesh area covered by particles",false,1.,"double",cmd);
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
    int tMax = iterationsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    int M = chainLengthArg.getValue(); 
    int nSources = sourcesArg.getValue();  
    int sampleNo = sampleNoArg.getValue();
    string meshName = meshSwitchArg.getValue();
    string savePath = savePathArg.getValue();
    double dt = deltaTArg.getValue();
    double areaFraction = areaFractionArg.getValue();
    double temperature = temperatureArg.getValue();
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = dangerousSwitch.getValue(); //not used right now

    shared_ptr<closedMeshSpace> meshSpace=make_shared<closedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
   
    double area = meshSpace->getArea();
    double sigma = 2*sqrt(areaFraction*area/(N*M_PI));
   
    if(programBranch >0)
        meshSpace->useSubmeshingRoutines(true,sigma,dangerous);

    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    configuration->setVerbose(verbose);
    configuration->setSpace(meshSpace);
    cout << "Area: " << meshSpace->getArea() << endl;

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,sigma);
    if(programBranch >= 1)
        configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly in a small space. Similarly, set random velocities in the tangent plane
    noiseSource noise(reproducible);
    configuration->setRandomParticlePositions(noise);
    configuration->setMaxwellBoltzmannVelocities(noise,temperature);

    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,sigma);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);

    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);

    //for now we fix timescale (tau) at 1.0
    double tau = 1.0;
    cout << "number of N-H masses: " << M << endl; 
    shared_ptr<noseHooverNVT> NVTUpdater=make_shared<noseHooverNVT>(dt, temperature, tau, M);
    NVTUpdater->setModel(configuration);
    simulator->addUpdater(NVTUpdater,configuration);
    
    profiler timer("various parts of the code");

    //by default, the simpleModelDatabase will save euclidean positions, mesh positions (barycentric + faceIdx), and particle velocities. See constructor for saving forces and/or particle types as well
    string meshFileName= meshName.substr(meshName.find_last_of("/") + 1);
    string extension = ".off";

    meshFileName = meshFileName.substr(0, meshFileName.size() - extension.size());
    
    char outputFileName[512];
    sprintf(outputFileName, "%s/NVT_N%i_areaFraction%.3f_sigma%.3f_tMax%i_mesh_%s_T%.2f_sample%i.h5",savePath.c_str(),N,areaFraction,sigma,tMax, meshFileName.c_str(),temperature,sampleNo);
    std::cout << "Creating new HDF5 file: " << outputFileName << std::endl;
    shared_ptr<simpleModelDatabase> saveState;
    saveState = make_shared<simpleModelDatabase>(N,outputFileName,fileMode::replace); 
    saveState->writeState(configuration,0.0);

    cout << "starting temp: " << NVTUpdater->getTemperatureFromKE() << endl; 
  
    ofstream temperatureFile("nvtTemperatures.csv"); 

    double density = N/area; 
    double rhoSigma2 = density*sigma*sigma; 
    
    double sumBPoR = 0;
    int counter = 1;
   
    string distancesFilename = "../sphere_data/NVT_N" + to_string(N) + "_a"+to_string(sigma) + "_distances.csv";
    std::ofstream distancesFile(distancesFilename);
    double largeDouble = 10000.0;
    int bonusTime = 2; 

    for (int ii = 0; ii < bonusTime*tMax; ++ii)
        {
        timer.start();
        simulator->performTimestep();
        timer.end();
        if(ii%saveFrequency == saveFrequency-1)
            {
            saveState->writeState(configuration,dt*ii);
            double nowTemp = NVTUpdater->getTemperatureFromKE();
            printf("step %i T %f \n",ii,nowTemp);
            temperatureFile << ii << ", " << nowTemp << "\n";
	    }
	if ((ii > tMax) && (ii%saveFrequency == saveFrequency-1))
       	    {
            double nowTemp = NVTUpdater->getTemperatureFromKE();
	    vector<double> stress(9, 0.0); 
	    simulator->computeMonodisperseStress(stress);
	    //stress trace is pressure, kb = 1, get temperature from average KE, and density is number density above
            sumBPoR += (stress[0]+stress[4]+stress[8])/(density*nowTemp);
	    counter+=1; 
            
	    //collect a set of distances every so often for later comparison
	    vector<meshPosition> positions = configuration->positions;
            //the below is necessary so that we don't accidentally cut off distances
            meshSpace->setNewSubmeshCutoff(largeDouble);
            meshSpace->useSubmeshingRoutines(false); 
            //use one particle as reference so number based normalization makes sense later.  
            int runningSum = 0;  
            
	    bool firstFlag = true; 
    	    //below should always get every distance available	    
	    for (int kk = 0; kk < nSources-1; kk ++) 
	        {
		vector<double> tempDistList;
                vector<vector3> startTangents;
                vector<vector3> endTangents;
		//don't want to count self distance (+1) or any previously calculated distances (+kk)	
	        vector<meshPosition>::const_iterator first = positions.begin() + 1 + kk;
                vector<meshPosition>::const_iterator last = positions.end(); 
                vector<meshPosition> targets(first, last);
                tempDistList.reserve(positions.size()); //largely just readability, a resize call happens within distance 
	        meshSpace->distance(positions[kk], targets, tempDistList, startTangents, endTangents);
	        //as usual, we only care about distances, so we could just discard startTangents and endTangents
		runningSum += tempDistList.size();
	        for (int jj = 0; jj < tempDistList.size(); jj ++)
       	            { 
		    if (firstFlag) 
		        {
                        distancesFile << tempDistList[jj];
			firstFlag = false;
		        }
		    else distancesFile << ", " << tempDistList[jj]; 
		    }
	        }
            distancesFile << "\n";
	    meshSpace->useSubmeshingRoutines(true, sigma, dangerous);
	    }
        }
    
    distancesFile.close(); 
    temperatureFile.close(); 
    
    timer.print();

   

    return 0;
    };
