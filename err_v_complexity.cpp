#include "std_include.h"
#include <tclap/CmdLine.h>
#include <iostream>
#include <fstream>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "cellListNeighborStructure.h"
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

using namespace TCLAP;
int main(int argc, char*argv[])
    {
    std::cout << "parsing command line." << std::endl;
    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("Setting up time vs. number of particles test...",' ',"V0.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> particlesArg("n", "N", "number of particles", false, 1000,"int", cmd);
    ValueArg<int> numDistsArg("i", "testNumber", "number of random distances to average over", false, 100, "int", cmd);  
    ValueArg<string> meshSwitchArg("m","meshSwitch","mesh goal surface id",false,"sphere","string",cmd);
    ValueArg<bool> verboseArg("v", "verbose", "verbosity", false, true, "bool", cmd); 
    SwitchArg submeshSwitch("z", "submeshed", "whether or not to use submesh", cmd, false); 
    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
     
    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int N = particlesArg.getValue();
    int numDists = numDistsArg.getValue();  
    string sName = meshSwitchArg.getValue();
    bool submeshed = submeshSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();

    bool verbose = true; // just always be verbose for tests 

    //for testing, particles initialized randomly. redo at the start of each check 
    //requires restarting the configuration 
    noiseSource noise(reproducible);
    
    //indicate submeshing or not     
    //create key simulation primitives 
    shared_ptr<simulation> simulator = make_shared<simulation>();
    
    string filename = "./"+sName+"_error_v_complexity.csv"; 
    std::ofstream complexity_tdf(filename); //tdf = tangents/distances/forces

    int nSurfaces = 13; 
    double largeDisplacement = 0.5;
    double smallDisplacement = 0.05;
    vector<vector3> emptyTransport = {};

    for (int s = 1; s < nSurfaces; s++) 
        {
	string meshName = "../err_v_c_meshes/"+sName + "_mesh_" +to_string(s)+".off"; 
	cout << "Checking remesh " << s << "." << endl;
        shared_ptr<triangulatedMeshSpace> meshSpace = make_shared<triangulatedMeshSpace>(); 
        meshSpace->loadMeshFromFile(meshName,verbose);
        double area = totalArea(meshSpace->surface);
        double particleSize = 2*sqrt(area/(N*M_PI)); //area fraction fixed at 1
        meshSpace->useSubmeshingRoutines(submeshed, particleSize);      
        int nv = meshSpace->surface.number_of_vertices();

	shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,particleSize);
        simulator->addForce(pairwiseForce);

	shared_ptr<simpleModel> configuration = make_shared<simpleModel>(N);
        configuration->setSpace(meshSpace);
        configuration->setRandomParticlePositions(noise);
        
        pairwiseForce->setModel(configuration);         

        simulator->setConfiguration(configuration);
	
	int numberTargets = 10; //breaks up the calculation to use a repeated source for numberTargets distances, saving time
			        //by re-using the sequence tree -- all-to-all is super expensive. We still get the same number of samples
				//if numDists%numTargets = 0.

	for (int jj = 0; jj < numDists/numberTargets; ++jj)
                {
                int index1 = noise.getInt(0,N-1);
		int i2; 
		vector<int> indices2;	
		while ((indices2.size() < numberTargets) && (N>1)) 
		    {
		    i2 = noise.getInt(0,N-1);
		    if (i2 != index1) indices2.push_back(i2);
                    }
		meshPosition source = configuration->positions[index1];
		vector<meshPosition> targets = {};
		for (int i: indices2) targets.push_back(configuration->positions[i]);
                
	        vector<vector3> displacementDirections;
	        //now combine mesh position weights of target points & face index of source point to get heading direction
		for (meshPosition mp: targets) 
	            {
		    meshPosition dispDummy; 
		    //because it's convenient, using the effectively random weights of the target mesh positions
		    //as ways to get random headings out from the source
		    dispDummy.x = mp.x;
		    dispDummy.faceIndex = source.faceIndex;
                    
		    vector<meshPosition> convertThese = {source, dispDummy}; 
		    
		    double3 emptyD3;
		    emptyD3.x = 0;
		    emptyD3.y = 0;
		    emptyD3.z = 0; 

                    vector<double3> converted = {emptyD3, emptyD3};
		    meshSpace->meshPositionToEuclideanLocation(convertThese, converted);
		    point3 sPoint(converted[0].x, converted[0].y, converted[0].z); 
		    point3 dPoint(converted[1].x, converted[1].y, converted[1].z);
		    displacementDirections.push_back(normalize(vector3(sPoint, dPoint))); 
		    }	    

		vector<double> distances;
                vector<vector3> startPath;
                vector<vector3> endPath;
                vector<vector3> displacements = {};
		meshSpace->distance(source,targets,distances,startPath,endPath);
                for (int ii = 0; ii < numberTargets; ii++) displacements.push_back(pairwiseForce->pairwiseForce(startPath[ii], distances[ii]));
                configuration->fillEuclideanLocations();
		for (int ii = 0; ii < numberTargets; ii++)
		    {	
		    vector3 pt = startPath[ii];
		    double3 p1 = configuration->euclideanLocations[index1]; 
		    double3 p2 = configuration->euclideanLocations[indices2[ii]];
		    complexity_tdf << nv << ", ";  
		    complexity_tdf << p1.x << ", " << p1.y << ", " << p1.z << ", ";
		    complexity_tdf << p2.x << ", " << p2.y << ", " << p2.z << ", ";
                    complexity_tdf << distances[ii] << ", "; 
		    complexity_tdf << pt.x() << ", " << pt.y() << ", " << pt.z() << ", ";
		    
		    vector3 displacement = displacementDirections[ii];
		    vector<meshPosition> drMeshPositions;
		    vector<double3> drs;
		    double3 dummyd3;
		    dummyd3.x = 0; 
		    dummyd3.y = 0;
		    dummyd3.z = 0;
		    drs.push_back(dummyd3); 
		    drs.push_back(dummyd3); 

		    meshPosition sourceCopy1;
		    meshPosition sourceCopy2; 
		    sourceCopy1.x = source.x;
		    sourceCopy1.faceIndex = source.faceIndex;
		    sourceCopy2.x = source.x;
		    sourceCopy2.faceIndex = source.faceIndex;
		    vector3 dLarge = largeDisplacement*displacement;
		    vector3 dSmall = smallDisplacement*displacement;
                    
		    meshSpace->transportParticleAndVectors(sourceCopy1, dLarge, emptyTransport);
                    meshSpace->transportParticleAndVectors(sourceCopy2, dSmall, emptyTransport);
	            drMeshPositions.push_back(sourceCopy1); 
		    drMeshPositions.push_back(sourceCopy2); 
                   
		    meshSpace->meshPositionToEuclideanLocation(drMeshPositions, drs); 
                   
		    double3 ldr = drs[0];
		    double3 sdr = drs[1]; 
		    complexity_tdf << displacement.x() << ", " << displacement.y() << ", " << displacement.z() << ", "; 
		    complexity_tdf << ldr.x << ", " << ldr.y << ", " << ldr.z << ", ";
		    complexity_tdf << sdr.x << ", " << sdr.y << ", " << sdr.z << "\n"; 
		     
		    }
		}

        }

    };
