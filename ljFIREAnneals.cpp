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
#include "lennardJones.h"
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

void saveNeighbors(vector<vector<int>> neighbors, string neighborFileName, int N)
    {
    ofstream neighborsFile;
    neighborsFile.open(neighborFileName, std::ios_base::app);
    for (int i = 0; i < N; i++) 
        {
        if (neighbors[i].size() > 0) 
	    {
	    neighborsFile << neighbors[i][0]; 
            }
	for (int j = 1; j < neighbors[i].size(); j++)
	    { 
            neighborsFile << ", " << neighbors[i][j]; 
    	    }
	neighborsFile << "\n";
	}
     }


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
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    
    ValueArg<double> afArg("f","areaFraction","packing fraction if every particle takes up circular space",false,1.,"double",cmd);
    ValueArg<double> omegaArg("w", "helicity", "helicity of mesh surface", false, 1.0, "double", cmd);  
    ValueArg<double> timestepArg("e", "timestepSize", "size of a timestep (not applicable for FIRE)", false, 0.1, "double", cmd);
    ValueArg<double> initializationRangeArg("a", "initializationRange", "range around origin to initialize in", false, 1.0, "double", cmd);

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);
    SwitchArg tangentialSwitch("c", "tangentialBCs", "use tangential BCs for open meshes", cmd, false); 
    SwitchArg trajectorySaveSwitch("p", "writeTrajectory", "write trajectory to file", cmd, false); 
    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int N = particleNumberSwitchArg.getValue();
    int totalAnneals = numberAnnealsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double areaFraction = afArg.getValue();
    double omega = omegaArg.getValue();
    double dt = timestepArg.getValue();
    double initializationRange = initializationRangeArg.getValue(); 
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool tangentialBCs = tangentialSwitch.getValue();
    bool saveTrajectory = trajectorySaveSwitch.getValue(); 
 
    profiler timer("various parts of the code");
    
    string inputMeshName = meshName.substr(14,15);
    string trajectoryFilename = "../bulk_silo_data/ljFIRE_"+to_string(N)+"_trajectories/" + inputMeshName + "_N_"+to_string(N) + "_trajectory.nc"; 
    cout << trajectoryFilename <<endl;
   
    string neighborFilename = "../bulk_silo_data/ljFIRE_"+to_string(N)+"_trajectories/neighbors_"+to_string(N)+"_"+to_string(areaFraction)+".csv"; 

    simpleModelDatabase saveState(N,trajectoryFilename,NcFile::Replace);
    ofstream neighborFile(neighborFilename);
       
    cout << "saving? " << saveTrajectory << endl;
    double fNorm, fMax, energyState; 
    int step = 1;//purely for printouts
    cout << "starting annealing loop" << endl;
    int GDSteps = 2000;   
    int FIRESteps = 500;     
    double cutoffSigma = 2.5;

    double maximumInteractionRange = 1; //placeholder
    shared_ptr<triangulatedMeshSpace> meshSpace = make_shared<triangulatedMeshSpace>(); //placeholder
    
    //!!!
    meshSpace->loadMeshFromFile(meshName,true);
    //!!!

    //for testing, just initialize particles randomly in a small space. Similarly, set random velocities in the tangent plane
    noiseSource noise(reproducible);

    vector<meshPosition> finalPositions;  

    for (int annealStep = 0; annealStep < totalAnneals; annealStep++) 
	{
	cout << "\nnew iteration, "<< annealStep << " updating..."<< endl;	
        meshSpace->updateMeshSpanAndTree(false);
	meshSpace->useTangentialBCs = tangentialBCs;
	cout << "setting positions are euclidean to false (?)" << endl;
        meshSpace->positionsAreEuclidean = false;

        double area = totalArea(meshSpace->surface);
        maximumInteractionRange = 2*sqrt(areaFraction*area/(N*M_PI));

	meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,false);

	shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
        configuration->setVerbose(verbose);
        configuration->setSpaceAndTMeshSpace(meshSpace);
 
        shared_ptr<simulation> simulator=make_shared<simulation>();
        simulator->setConfiguration(configuration);

	shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);     
        shared_ptr<fireMinimization> fireEnergyMinimizer=make_shared<fireMinimization>();
	//becuase we're using FIRE with LJ, it's critical that the max timestep size is at most very small! here it's set to dt/10. 
        fireEnergyMinimizer->setFIREParameters(1, dt, 0.99, dt/10, 1e-10, 1.1, 0.9, 0.9, 4, 1e-10,0.00);
	energyMinimizer->setModel(configuration);
        fireEnergyMinimizer->setModel(configuration); 
	
	shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.00,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
        shared_ptr<lennardJones> LJForce = make_shared<lennardJones>(1.00, maximumInteractionRange, true); 
        LJForce->setModel(configuration);
        LJForce->cutoffCoefficient = cutoffSigma;
	pairwiseForce->setModel(configuration);
    
        shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);

	configuration->setNeighborStructure(cellList); 
	configuration->setRandomMeshPositionsNearZero(noise, initializationRange);
        configuration->setMaxwellBoltzmannVelocities(noise, 0);

	if (saveTrajectory) saveState.writeState(configuration,0.0);
		
	simulator->addUpdater(energyMinimizer, configuration);       
        simulator->addForce(pairwiseForce);
        
   	
	printf("GD Cooling\n"); 
	for (int ii = 0; ii < GDSteps; ii++)
            {
	    timer.start();
            simulator->performTimestep();
            timer.end();	    
            
	    energyState = pairwiseForce->computeEnergy();
	    fNorm = energyMinimizer->getForceNorm();

            if(step%saveFrequency == saveFrequency-1)
                {
                if (saveTrajectory) saveState.writeState(configuration,step); 
                fMax = energyMinimizer->getMaxForce();

	        configuration->findNeighbors(maximumInteractionRange); 
	        saveNeighbors(configuration->neighbors, neighborFilename, N);	
		
		printf("omega %.4g step %i fN %.4g fM %.4g E %.4g\n",omega, step,fNorm,fMax, energyState);
		}

	    step++;
            };
        //reset updater/force computers for LJ force, using FIRE to relax quickly
	simulator->clearUpdaters();
        simulator->clearForceComputers();        
	
	simulator->addUpdater(fireEnergyMinimizer,configuration);
	meshSpace->setNewSubmeshCutoff(cutoffSigma*maximumInteractionRange);

	simulator->addForce(LJForce);
	simulator->computeForces(); 
 
        /*	
        cout << "printing first lj forces" << endl;
        for (auto f: configuration->forces) cout << f << "\n"; 
        */

	printf("FIRE cooling\n"); 
	for (int ii = 0; ii < FIRESteps; ii++) 
	    {
	    timer.start();
	    simulator->performTimestep();
	    timer.end();
	    if(step%saveFrequency == saveFrequency-1)
                {
                if (saveTrajectory) saveState.writeState(configuration,step);
                fMax = fireEnergyMinimizer->getMaxForce();
                fNorm = fireEnergyMinimizer->getForceNorm();
		energyState = LJForce->computeEnergy();
		
		configuration->findNeighbors(cutoffSigma*maximumInteractionRange);
		saveNeighbors(configuration->neighbors, neighborFilename, N);
		printf("omega %.4g step %i fN %.16g fM %.16g E %.16g\n",omega, step,fNorm,fMax, energyState);
                }
            step++;     
	    }    

        string minimumFilename = "../bulk_silo_data/ljFIRE_"+to_string(N)+"_minimals/" + inputMeshName + "_" + to_string(N) + "_minimal_anneal"+to_string(annealStep)+".nc";
	
	//save the configuration & its voronoi positions at the end of the relaxation	    
        simpleModelDatabase minState(N, minimumFilename, NcFile::Replace); //don't use posToSave.size(), it's three times as large as it should be 
        
	minState.writeState(configuration,step); 
	configuration->fillEuclideanLocations();
        string mposFilename = "../bulk_silo_data/ljFIRE_"+to_string(N)+"_voronois/" + inputMeshName + to_string(N) + "_voronoi_positions"+to_string(annealStep)+".csv";
	writeVoronoiPositions(configuration->positions, configuration->euclideanLocations, mposFilename, N);  
        if (annealStep == totalAnneals-1) 
	    {
	    cout << "writing final positions for r/d comparison" << endl;
            finalPositions = configuration->positions;
            }
	};

    timer.print();    

    point3 origin(0,0,0);
    pmpFaceLocation originBary = PMP::locate(origin,meshSpace->surface);
    meshPosition originMeshPos;
    originMeshPos.x = point3(originBary.second[0],originBary.second[1],originBary.second[2]);
    originMeshPos.faceIndex = int(originBary.first);
    vector<double> allDists;
    vector<vector3> startTangents;
    vector<vector3> endTangents;

    //the below is necessary so that we don't accidentally cut off distances, since we actually want the genuine 
    //maximum
    double largeDouble = 10000.0;
    meshSpace->setNewSubmeshCutoff(largeDouble);
    meshSpace->useSubmeshingRoutines(false);

    meshSpace->distance(originMeshPos, finalPositions, allDists, startTangents,endTangents);
    sort(allDists.begin(), allDists.end(), greater<double>());
    int nToSum = 5;
    double sum=0;
    for (int i = 0; i < nToSum; i++) 
        {
        sum+=allDists[i];
	cout << allDists[i] << endl;
        }
    double R = sum/nToSum; // we could average instead from e.g. 5 maximum distances or n/50 maximum distances

    double omegaR = omega*R;
    double RoverD = R/maximumInteractionRange;
    cout << "R: " << R << ", d: " << maximumInteractionRange << endl; 
    cout << "omegaR = " << omegaR << endl;
    cout << "R over d = " << RoverD << endl;
     
    ofstream plotPointsFile;
    plotPointsFile.open("../bulk_silo_data/ljFIRE_omegaR_Rdivd" + to_string(N) + ".csv", ios::app);
    if (!plotPointsFile)
        {
        plotPointsFile = ofstream("../bulk_silo_data/ljFIRE_omegaR_Rdivd.csv");
        plotPointsFile << "Omega, R, OmegaR, RoverD" << endl;	
        }
    plotPointsFile << omega << ", " << R <<", "<< omegaR << ", " << RoverD << ", " << "\n";  
    plotPointsFile.close();  

    return 0;
    };

