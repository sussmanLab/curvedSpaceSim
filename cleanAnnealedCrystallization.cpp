#include "std_include.h"
#include <tclap/CmdLine.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "fireMinimization.h"
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
    bool outFlag = false;
    int counter = 0;
    double bound = 2.0; 
    for (meshPosition p: positions)
        {
        double w1 = p.x.x();
	double w2 = p.x.y(); 
	double w3 = p.x.z(); 
        //it's funny to look for the right structure for "not greater than,"
        //as people are quick to quip that you can use less than (ignoring nan case)
        if (!(w1 < bound && w1 > -bound)) outFlag = true;
        if (!(w2 < bound && w2 > -bound)) outFlag = true;
	if (!(w3 < bound && w3 > -bound)) outFlag = true;
	if (outFlag) 
	    {
	    cout << "Step " << step << endl;  
	    cout << "particle " << counter << endl;
	    cout << "particle weights: "; 
	    printPoint(p.x); 
	    printf("\n");	   
            ERRORERROR("Bad particle found.");
	    }
	counter++;
	}
    } 

void checkEucOutOfBounds(vector<double3> positions, vector<meshPosition> baryCoords, int step)
    {
    bool outFlag = false;
    int counter = 0;
    double bound = 1000.0;
    for (double3 p: positions)
        {
        double x = p.x;
        double y = p.y;
        double z = p.z;
        //it's funny to look for the right structure for "not greater than,"
        //as people are quick to quip that you can use less than (ignoring nan case)
        if (!(x < bound && x > -bound)) outFlag = true;
        if (!(y < bound && y > -bound)) outFlag = true;
        if (!(z < bound && z > -bound)) outFlag = true;
        counter++;
	if (outFlag) break; 
        }
    if (outFlag)
        {
        cout << "Step " << step << endl;
        cout << "broken particle " << counter << endl; 
	cout << "particle positions: \n";
	int c = 0; 
        for (double3 p: positions)
	    {
            printf("%i ", c);
	    printf("{%f,%f,%f}",p.x,p.y,p.z);
	    printf("\n");
	    c++;
	    }
	
	c=0;
	cout << "barycentric coords: \n"; 
        for (meshPosition b: baryCoords)
            {
            printf("%i ", c);
            printPoint(b.x);
            printf("\n");
	    c++; 
            }

        ERRORERROR("Bad particle Euclidean location found.");
        }
    }

vector<int> findExclusions(vector<meshPosition> &positions, triangleMesh &theMesh)
    {
    vector<int> exclusions;
    exclusions.reserve(positions.size()/2); //way more than it needs
    for (int ii = 0; ii < positions.size(); ii++)
        {
        faceIndex iiFace = faceIndex(positions[ii].faceIndex);
        //is_border only works with halfedges, so we need to dig out the associated edges
        halfedgeIndex hf(theMesh.halfedge(iiFace));
        for(halfedgeIndex hi : halfedges_around_face(hf, theMesh))
            {
            edgeIndex ei = theMesh.edge(hi);
            if(theMesh.is_border(ei))
                {
                exclusions.push_back(ii);
                break;
                }
            }
        }
    return exclusions;
    }

void writeVoronoiPositions(vector<meshPosition> meshPositions, vector<double3> R3Positions, string filename, int N)
    {
    std::ofstream voronoiPosFile(filename);

    for (int ii = 0; ii < N; ++ii)
        {
        int fIndex = meshPositions[ii].faceIndex;
        double3 p = R3Positions[ii];
        voronoiPosFile << fIndex << ", " << p.x << ", " << p.y << ", " << p.z << "\n";
        }
    };


using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in curved space!",' ',"V0.9");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> saveFrequencyArg("s","saveFrequency","how often a file gets updated",false,100,"int",cmd); 
    ValueArg<int> numberAnnealsArg("i","nAnneals","number of anneals to do",false, 10, "int",cmd);
    ValueArg<int> chainLengthArg("l", "M", "length of N-H chain", false, 2, "int", cmd);
    ValueArg<int> programBranchArg("z", "branch", "program branch", false, 0, "int", cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<double> densityArg("d","density","packing fraction if every particle takes up circular space",false,1.,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,1.0,"double",cmd);
    ValueArg<double> stiffnessArg("q", "stiffness", "harmonic stiffness", false, 1.0, "double", cmd); 
    ValueArg<double> omegaArg("w", "helicity", "helicity of mesh surface", false, 1.0, "double", cmd);  
    ValueArg<double> timestepArg("e", "timestepSize", "size of a timestep (not applicable for FIRE)", false, 0.1, "double", cmd);

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);
    SwitchArg tangentialSwitch("c", "tangentialBCs", "use tangential BCs for open meshes", cmd, false); 

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int N = particleNumberSwitchArg.getValue();
    int M = chainLengthArg.getValue();
    int totalAnneals = numberAnnealsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    int programBranch = programBranchArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double density = densityArg.getValue();
    double temperature = temperatureArg.getValue();
    double stiffness = stiffnessArg.getValue();  
    double omega = omegaArg.getValue();
    double dt = timestepArg.getValue(); 
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool tangentialBCs = tangentialSwitch.getValue(); 


    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);

    double area = totalArea(meshSpace->surface);
    double maximumInteractionRange = 2*sqrt(density*area/(N*M_PI));

    meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,false);
    meshSpace->useTangentialBCs = tangentialBCs;

    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    configuration->setVerbose(verbose);
    configuration->setSpaceAndTMeshSpace(meshSpace);

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
    configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly in a small space. Similarly, set random velocities in the tangent plane
    noiseSource noise(reproducible);
   
    /* 
    if (positionsFile != "none") {        
       configuration->setMeshPositionsFromR3File(positionsFile, meshSpace->surface);
    }*/
    configuration->setRandomParticlePositions(noise);

 
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

    cout << "creating input mesh name" << endl;    
    string inputMeshName = meshName.substr(14,15);
    cout << "input mesh name " << inputMeshName << endl; 
    string trajectoryFilename = "../bulk_silo_data/siloMinimization_" + to_string(N) + "_" + inputMeshName + ".nc"; 
    simpleModelDatabase saveState(N,trajectoryFilename,NcFile::Replace);
    saveState.writeState(configuration,0.0);

    double fNorm, fMax, energyState; 
    fNorm=1;
    int heatingLength = 1000; 
    double runningMinimum = 2*(pairwiseForce->computeEnergy()); 
    string minimumFilename = "../bulk_silo_data/minimals/siloMinimization_" + to_string(N) + "_" + inputMeshName + "_minconfig.nc"; 
    int step = 1;//purely for printouts
    double forceNormCutoff = 1e-9;
    cout << "starting annealing loop" << endl;
    vector<double3> minR3Positions;
    vector<meshPosition> minMeshPositions;  
    int maxDescentSteps = 100000; //number of fire steps to take at most, each 100 long  
    
    //ofstream forceFile("forces_"+to_string(N)+"_"+inputMeshName+".csv");

    for (int annealStep = 0; annealStep < totalAnneals; annealStep++) 
	{
	shared_ptr<fireMinimization> fireEnergyMinimizer=make_shared<fireMinimization>();
        fireEnergyMinimizer->useFWithExclusions=true;
	shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
	
	simulator->clearUpdaters();
	vector<int> preExclusions = findExclusions(configuration->positions, meshSpace->surface);	
	if (programBranch ==1 )
       	    {
	    printf("fire minimization\n");  
	    fireEnergyMinimizer->setModel(configuration); 
   	    simulator->addUpdater(fireEnergyMinimizer,configuration);
            fNorm = fireEnergyMinimizer->getForceNormWithExclusions(preExclusions);
	    }
        else 
            {
	    printf("gradient descent minimization\n");
	    energyMinimizer->setModel(configuration);
	    simulator->addUpdater(energyMinimizer, configuration);
            fNorm = energyMinimizer->getForceNormWithExclusions(preExclusions);
	    }

	for (int ii = 0; ii < maxDescentSteps; ii++)
            {
	
	    checkOutOfBounds(configuration->positions,step);
	    configuration->fillEuclideanLocations(); 
            checkEucOutOfBounds(configuration->euclideanLocations, configuration->positions, step);	
	    timer.start();
            simulator->performTimestep();
            timer.end();	    
            energyState = pairwiseForce->computeEnergy();
	    
	    vector<int> exclusions = findExclusions(configuration->positions, meshSpace->surface);
            
	    if (programBranch == 1) fNorm = fireEnergyMinimizer->getForceNormWithExclusions(exclusions);
	    else fNorm = energyMinimizer->getForceNormWithExclusions(exclusions);
	    
	    if (fNorm <= forceNormCutoff) 
	        {
		printf("Reached force norm threshold.\n"); 
		if (programBranch == 1) fMax = fireEnergyMinimizer->getMaxForceWithExclusions(exclusions);
                else fMax = energyMinimizer->getMaxForceWithExclusions(exclusions);
		printf("step %i fN %.4g fM %.4g E %.4g\n",step,fNorm,fMax, energyState);
		saveState.writeState(configuration,step); 
		//forceFile << fNorm << ", " << fMax << endl;
		break; 
		}

            if(step%saveFrequency == saveFrequency-1)
                {
                saveState.writeState(configuration,step); 
                if (programBranch == 1) fMax = fireEnergyMinimizer->getMaxForceWithExclusions(exclusions);
                else fMax = energyMinimizer->getMaxForceWithExclusions(exclusions);
		printf("step %i fN %.4g fM %.4g E %.4g\n",step,fNorm,fMax, energyState);
                //forceFile << fNorm << ", " << fMax << endl;
		}
	    step++;
	    if (ii == maxDescentSteps-1) printf("Reached max cooling length; stopping."); 
            };
        
	//now, save the minimum energy configuration
	if (energyState < runningMinimum) 
	    {
	    //if this is our lowest energy configuration, save the energy value and the configuration itself	    
            simpleModelDatabase minState(N, minimumFilename, NcFile::Replace); //don't use posToSave.size(), it's three times as large as it should be 
            //because it's defined local to this conditional, above will overwrite existing files
	    runningMinimum=energyState;
	    minState.writeState(configuration,step); 
	    configuration->fillEuclideanLocations();
	    checkEucOutOfBounds(configuration->euclideanLocations, configuration->positions, step);
	    minR3Positions = configuration->euclideanLocations; 
	    minMeshPositions = configuration->positions;  
	    }
        
	shared_ptr<noseHooverNVT> NVTUpdater=make_shared<noseHooverNVT>(dt, temperature, 1.0, M);
        NVTUpdater->setModel(configuration); 
	simulator->clearUpdaters();
        
	configuration->setMaxwellBoltzmannVelocities(noise,temperature);
        simulator->addUpdater(NVTUpdater,configuration);
	//now, heat the system up to get a new configuration 
	for (int ii = 0; ii < heatingLength; ++ii)
	    {
	    vector<meshPosition> cPositions = configuration->positions;
	    vector<vector3> cVelocities = configuration->velocities; 
                        checkOutOfBounds(configuration->positions,step);
            configuration->fillEuclideanLocations();
            checkEucOutOfBounds(configuration->euclideanLocations, configuration->positions, step);
	    simulator->performTimestep();
	    step++;
	    if(step%saveFrequency == saveFrequency-1)
                {
                saveState.writeState(configuration,step);
		double nowTemp = NVTUpdater->getTemperatureFromKE();
                printf("step %i T %f \n",step,nowTemp);
                }
	    }
	};
    timer.print();    
    
    
    //substr size is meant to include the omega section of silo_meshes mesh names 
    string mposFilename = "../bulk_silo_data/voronois/silo_" + to_string(N) + inputMeshName + "_voronoi_positions_clean.csv";
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
    cout << "Minimal energy found: " << runningMinimum << endl;
    cout << "R: " << R << ", d: " << maximumInteractionRange << endl; 
    cout << "omegaR = " << omegaR << endl;
    cout << "R over d = " << RoverD << endl;
     
    ofstream plotPointsFile;
    plotPointsFile.open("../bulk_silo_data/clean_omegaR_Rdivd" + to_string(N) + ".csv", ios::app);
    if (!plotPointsFile)
        {
        plotPointsFile = ofstream("../bulk_silo_data/clean_omegaR_Rdivd.csv"); 
        }
    plotPointsFile << omega << ", " <<omegaR << ", " << RoverD << ", " << minimumFilename << "\n";  
    plotPointsFile.close();  

    return 0;
    };

