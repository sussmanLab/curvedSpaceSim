#include "std_include.h"
#include <tclap/CmdLine.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "simpleModel.h"
#include "noseHooverNVT.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "simpleModelDatabase.h"
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

void checkOutOfBounds(vector<meshPosition> positions, int step) 
    {
    bool nanFlag = false;
    for (meshPosition p: positions)
        {
        double w1 = p.x.x();
	double w2 = p.x.y(); 
	double w3 = p.x.z(); 
        //it's funny to look for the right structure for "not greater than,"
        //as people are quick to quip that you can use less than (ignoring nan case)
        if (!((w1 < 2) && (w1 > -1))) {nanFlag = true; break;}
        if (!((w2 < 2) && (w2 > -1))) {nanFlag = true; break;}
	if (!((w3 < 2) && (w3 > -1))) {nanFlag = true; break;}
	}
        if (nanFlag) cout << "Step " << step << endl;
        if (nanFlag) ERRORERROR("Bad particle found." );
    } 

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
    ValueArg<int> descentStepsArg("i","GDsteps","steps in each gradient descent",false, 10000, "int",cmd);
    ValueArg<int> chainLengthArg("l", "M", "length of N-H chain", false, 2, "int", cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<string> positionsFileArg("f", "positionsFile", "filename for xyz positions you want to load", false, "none", "string", cmd); 
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("e","dt","timestep size",false,.01,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,1.0,"double",cmd);
    ValueArg<double> stiffnessArg("q", "stiffness", "harmonic stiffness", false, 1.0, "double", cmd); 
    ValueArg<double> omegaArg("w", "helicity", "helicity of mesh surface", false, 1.0, "double", cmd);  

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);
    SwitchArg tangentialSwitch("c", "tangentialBCs", "use tangential BCs for open meshes", cmd, false); 

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int M = chainLengthArg.getValue();
    int descentSteps = descentStepsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    string meshName = meshSwitchArg.getValue();
    string positionsFile = positionsFileArg.getValue();
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    double temperature = temperatureArg.getValue();
    double stiffness = stiffnessArg.getValue();  
    double omega = omegaArg.getValue(); 
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool tangentialBCs = tangentialSwitch.getValue(); 


    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
    if(programBranch >0)
        meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,false);
    meshSpace->useTangentialBCs = tangentialBCs;

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

 
    //riesz potential order is scale, stiffness exponent
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(stiffness,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    //shared_ptr<rieszPotential> pairwiseForce = make_shared<rieszPotential>(energyScale,3.0);
    pairwiseForce->setModel(configuration);

    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);
        
    profiler timer("various parts of the code");

    vector<double> posToSave;
    getFlatVectorOfPositions(configuration,posToSave);
    
    string inputMeshName = meshName.substr(15,14); 
    string trajectoryFilename = "./siloMinimization_" + to_string(N) + "_" + inputMeshName + ".nc"; 
    simpleModelDatabase saveState(N,trajectoryFilename,NcFile::Replace);
    saveState.writeState(configuration,0.0);

    double fNorm,fMax, energyState; 
    fNorm=1;
    int totalAnneals = 20;
    int heatingLength = 1000; 
    double runningMinimum = 2*(pairwiseForce->computeEnergy()); 
    string minimumFilename = "../silo_data/siloMinimization_" + to_string(N) + "_" + inputMeshName + "_minconfig.nc";

    int step = 1;
    cout << "starting annealing loop" << endl;
    vector<double3> minR3Positions;
    vector<meshPosition> minMeshPositions;  
    meshPosition obviouslyWrongMeshPosition; 

    for (int annealStep = 0; annealStep < totalAnneals; annealStep++) 
	{
        shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);		
        energyMinimizer->setModel(configuration);
	simulator->clearUpdaters();
   	simulator->addUpdater(energyMinimizer,configuration);
	//gradient descent minimization	
        fNorm = energyMinimizer->getForceNorm();
        //before the cooling section of every heating/cooling cycle, we'll check in with the 
	//positions to see if anything has gone awry
	
        for (meshPosition pos: configuration->positions)
	{ cout << endl; printPoint(pos.x);}
        cout << endl;	

	for (int ii = 0; ii < descentSteps; ii++)
            {
	
	    checkOutOfBounds(configuration->positions,step);
	
	    timer.start();
            simulator->performTimestep();
            timer.end();	    

            if(step%saveFrequency == saveFrequency-1)
                {
                saveState.writeState(configuration,dt*step); 
		if(programBranch <2)
                    {
                    fNorm = energyMinimizer->getForceNorm();
                    fMax = energyMinimizer->getMaxForce();
		    energyState = pairwiseForce->computeEnergy();
                    printf("step %i fN %.4g fM %.4g E %.4g\n ",step,fNorm,fMax, energyState);
                    }
                else
                    printf("step %i \n",step);
               }
	    step++; 
            };
        
	//now, save the minimum energy configuration
        double currentMinEnergy = pairwiseForce->computeEnergy(); 
	if (currentMinEnergy < runningMinimum) 
	    {
	    //if this is our lowest energy configuration, save the energy value and the configuration itself	    
            simpleModelDatabase minState(posToSave.size(), minimumFilename, NcFile::Replace); 
            //because it's defined local to this conditional, above will overwrite existing files
	    runningMinimum=currentMinEnergy;
	    minState.writeState(configuration,dt*step); 
	    configuration->fillEuclideanLocations();
	    minR3Positions = configuration->euclideanLocations; 
	    minMeshPositions = configuration->positions;  
	    }
        
	shared_ptr<noseHooverNVT> NVTUpdater=make_shared<noseHooverNVT>(dt, temperature, 1.0, M);
        NVTUpdater->setModel(configuration); 
	simulator->clearUpdaters();
	cout << "setting mb velocities" << endl;
        configuration->setMaxwellBoltzmannVelocities(noise,temperature);
        simulator->addUpdater(NVTUpdater,configuration);
	//now, heat the system up to get a new configuration 
	for (int ii = 0; ii < heatingLength; ++ii)
	    {
	    checkOutOfBounds(configuration->positions,step);
            simulator->performTimestep();
	    step++;
	    if(step%saveFrequency == saveFrequency-1)
                {
                //this needs to be updated to new saveState style 
		//(see curvedSpaceNVTSim) but for now just proof of concept
                saveState.writeState(configuration,dt*step);
		double nowTemp = NVTUpdater->getTemperatureFromKE();
                printf("step %i T %f \n",step,nowTemp);
                }
	    }
	};
    
    
    timer.print();    
    
    
    //substr size is meant to include the omega section of silo_meshes mesh names 
    string mposFilename = "../silo_data/silo_" + to_string(N) + inputMeshName + "_voronoi_positions.csv";
    std::ofstream meshPosFile(mposFilename); 
    
    for (int ii = 0; ii < N; ++ii)
        {
	int fIndex = minMeshPositions[ii].faceIndex; 
        double3 p = minR3Positions[ii];
	meshPosFile << fIndex << ", " << p.x << ", " << p.y << ", " << p.z << "\n";
        }

    point3 origin(0,0,0);
    pmpFaceLocation originBary = PMP::locate(origin,meshSpace->surface);
    meshPosition originMeshPos;
    originMeshPos.x = point3(originBary.second[0],originBary.second[1],originBary.second[2]);
    originMeshPos.faceIndex = int(originBary.first);
    double maxDist = 0;
    vector<double> allDists;
    vector<vector3> startTangents;
    vector<vector3> endTangents;

    //the below is necessary so that we don't accidentally cut off distances, since we actually want the genuine 
    //maximum
    double largeDouble = 10000.0;
    meshSpace->setNewSubmeshCutoff(largeDouble);
    meshSpace->useSubmeshingRoutines(false);

    meshSpace->distance(originMeshPos, minMeshPositions, allDists, startTangents,endTangents);
    for (double dist: allDists)
        {
        maxDist = max(maxDist,dist);
        }
    double R = maxDist; // we could average instead from e.g. 5 maximum distances or n/50 maximum distances

    double omegaR = omega*R;
    double RoverD = R/maximumInteractionRange;
    cout << "R: " << R << ", d: " << maximumInteractionRange << endl; 
    cout << "omegaR = " << omegaR << endl;
    cout << "R over d = " << RoverD << endl;
   
    ofstream plotPointsFile;
    plotPointsFile.open("../silo_data/omegaR_Rdivd.csv", ios::app);
    if (!plotPointsFile)
        {
        plotPointsFile = ofstream("omegaR_Rdivd.csv"); 
        }
    plotPointsFile << omegaR << ", " << RoverD << ", " << minimumFilename << "\n";  
    plotPointsFile.close();  

    return 0;
    };
