#include "std_include.h"
#include <tclap/CmdLine.h>


#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "fireMinimization.h"
#include "velocityVerletNVE.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "lennardJones.h"
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
    ValueArg<int> recordArg("l", "recordNumber", "record number to pull from trajectory file", false, 0, "int", cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../bulk_silos/silo_omega0.004000_R2.000000_remesh.off","string",cmd);
    ValueArg<string> trajectoryFileArg("v", "trajectoryFile", "filename of the trajectory you want to load", false, "NULL", "string", cmd);
    ValueArg<double> initializationRangeArg("a","initializationRange","distance from origin to allow particles to initialize",false,1.,"double",cmd);
    ValueArg<double> deltaTArg("e","dt","timestep size",false,.01,"double",cmd);
    ValueArg<double> temperatureArg("t","T","temperature to set",false,.2,"double",cmd);
    ValueArg<double> afArg("f","areaFraction","area fraction for disks on the surface",false,1.0,"double",cmd);

    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg exclusionsSwitch("x", "fExclusions", "use force with exclusions", cmd, false);
    SwitchArg siloSwitch("o", "usingSiloMesh", "o for open boundary...", cmd, true); 
    
    //parse the arguments
    cmd.parse(argc, argv);
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int maximumIterations = iterationsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    int recordToLoad = recordArg.getValue();
    string trajectoryFile = trajectoryFileArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double initializationRange = initializationRangeArg.getValue();
    double dt = deltaTArg.getValue();
    double temperature = temperatureArg.getValue();
    double areaFraction = afArg.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool excludeBoundary = exclusionsSwitch.getValue();

    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,false);
    meshSpace->useSubmeshingRoutines(false);
    
    double area = totalArea(meshSpace->surface);
    double maximumInteractionRange = 2*sqrt(areaFraction*area/(N*M_PI));
    cout << "max interaction range " << maximumInteractionRange << endl;
   
    if(programBranch >=0) meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,false);
    meshSpace->useTangentialBCs = true;
 
    shared_ptr<simpleModel> configuration=make_shared<simpleModel>(N);
    configuration->setVerbose(false);
    configuration->setSpaceAndTMeshSpace(meshSpace);

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
    if(programBranch >= 0)
        configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly in a small space. Similarly, set random velocities in the tangent plane
    noiseSource noise(reproducible);
    //configuration->setRandomParticlePositions(noise);
    if (trajectoryFile == "NULL")
        {
	configuration->setRandomMeshPositionsNearZero(noise, initializationRange);
        configuration->setMaxwellBoltzmannVelocities(noise,temperature);
	}
    else 
        {
        simpleModelDatabase loadState(N,trajectoryFile,NcFile::ReadOnly);
        loadState.readState(configuration, recordToLoad);  
        }

    //shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,.5);
    shared_ptr<harmonicRepulsion> pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);
    double ljEnergyScale = 1.00;
    shared_ptr<lennardJones> LJForce = make_shared<lennardJones>(ljEnergyScale, maximumInteractionRange, true); 
    double cutoffSigma = 2.5;
    LJForce->cutoffCoefficient = cutoffSigma;

    shared_ptr<simulation> simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);

    shared_ptr<gradientDescent> energyMinimizer=make_shared<gradientDescent>(dt);
    shared_ptr<fireMinimization> energyMinimizerFire=make_shared<fireMinimization>();
    shared_ptr<velocityVerletNVE> nve=make_shared<velocityVerletNVE>(dt);
     
   //void setFIREParameters(int _maximumIterations,double _deltaT, double _alphaStart, double _deltaTMax, double _deltaTMin, double _deltaTInc, double _deltaTDec, double _alphaDec, int _nMin, double _forceCutoff, double _alphaMin = 0.75);
    double maxDT = dt/10; 
    energyMinimizerFire->setFIREParameters(1, dt, 0.99, maxDT, 1e-10, 1.1, 0.9, 0.9, 4, 1e-10,0.00); 
    energyMinimizerFire->useFWithExclusions = excludeBoundary;
    cout <<"continuing using FIRE minimization" << endl;
    energyMinimizerFire->setModel(configuration);
    simulator->addUpdater(energyMinimizerFire,configuration);
    
    profiler timer("various parts of the code");

    //by default, the simpleModelDatabase will save euclidean positions, mesh positions (barycentric + faceIdx), and particle velocities. See constructor for saving forces and/or particle types as well
    string minimizerName = "FIRE"; 
    simpleModelDatabase saveState(N,"./ljContinuation_" + to_string(N) +"_af" + minimizerName+to_string(areaFraction)+".nc",NcFile::Replace, true, false, true);
    saveState.writeState(configuration,0.0); 
    ofstream forceFile("./forces_"+to_string(N)+"_"+minimizerName+to_string(areaFraction)+".csv");
    
    ofstream neighborsFile("./continuationNeighbors_"+to_string(N)+minimizerName+to_string(areaFraction)+".csv"); 
    ofstream distanceFile("./continuationDistances_"+to_string(N)+minimizerName+to_string(areaFraction)+".csv");
    cout << "Using exclusions? " << excludeBoundary << endl;
    
    energyMinimizerFire->setModel(configuration);
    simulator->addUpdater(energyMinimizerFire,configuration);
    LJForce->setModel(configuration);
	    
    meshSpace->setNewSubmeshCutoff(cutoffSigma*maximumInteractionRange); 
    simulator->addForce(LJForce);
    simulator->computeForces();
    programBranch = 0;
    cout << "printing first lj forces" << endl;
    
    for (auto f: configuration->forces) cout << f << "\n"; 
    
    for (int ii = 0; ii < maximumIterations; ++ii)
        {	       
	timer.start();
	simulator->performTimestep();
        timer.end();
        
	double energy; 
	energy = LJForce->computeEnergy(); 
            
        //purely a neighbor saver, not related to main routine unless by accident	
	configuration->findNeighbors(cutoffSigma*maximumInteractionRange);
	vector<vector<int>> neighbors = configuration->neighbors;
	vector<vector<double>> neighborDistances = configuration->neighborDistances; 

        for (int i = 0; i < N; i++)
            {
	    if (neighbors[i].size() > 0) 
	        {
		neighborsFile << neighbors[i][0]; 
		distanceFile << neighborDistances[i][0];
		}
	    for (int j = 1; j < neighbors[i].size(); j++)
	        { 
                neighborsFile << ", " << neighbors[i][j]; 
                distanceFile << ", " << neighborDistances[i][j]; 
    	        }
	    neighborsFile << "\n";
            distanceFile << "\n";
	    }

        if(ii%saveFrequency == saveFrequency-1)
            {
            saveState.writeState(configuration,dt*ii);
	    vector<int> exclusions;
	    configuration->findBoundaryParticles(exclusions);
            if(programBranch == 1)
                {
                double fNorm,fMax;
		if (excludeBoundary)
		    {
                    fNorm = energyMinimizer->getForceNormWithExclusions(exclusions);
                    fMax = energyMinimizer->getMaxForceWithExclusions(exclusions);
		    }
		else
		    {
		    fNorm = energyMinimizer->getForceNorm();
                    fMax = energyMinimizer->getMaxForce();
                    }
		printf("step %i fN %g fM %g E %.16g\n",ii,fNorm,fMax, energy);
	        forceFile << ii << ", " << fNorm << ", " << fMax << "\n";	
		}
            else 
	        { 
                double fNorm,fMax;
		if (excludeBoundary)
                    {
                    fNorm = energyMinimizerFire->getForceNormWithExclusions(exclusions);
                    fMax = energyMinimizerFire->getMaxForceWithExclusions(exclusions);
                    }
                else
                    {
                    fNorm = energyMinimizerFire->getForceNorm();
                    fMax = energyMinimizerFire->getMaxForce();
                    }
                printf("step %i fN %g fM %g E %.16g\n",ii, fNorm,fMax, energy);
                forceFile << ii << ", " << fNorm << ", " << fMax << "\n";
                  
              }
            }
        }
    timer.print();

    return 0;
    };
