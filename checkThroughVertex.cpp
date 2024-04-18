#include "triangulatedMeshSpace.h"

int main(int argc, char*argv[]) 
    {
    string meshName = "../exampleMeshes/torus_isotropic_remesh.off"; 
    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,true);
	    
    //this pair goes directly through a vertex! which vertex? does it matter? it's a test!

    point3 weights1(0.4,0.3,0.3);
    point3 weights2(0.7,0.3,0);
    meshPosition source(weights1, 1);
    meshPosition target(weights2, 1);
    
    std::vector<double3> euc1s;
    std::vector<double3> euc2s;

    std::vector<meshPosition> sources;
    std::vector<meshPosition> targets;
    sources.push_back(source);
    targets.push_back(target);
    //nab r3 dispalcement by converting to euclidean 
    meshSpace->meshPositionToEuclideanLocation(sources, euc1s);
    meshSpace->meshPositionToEuclideanLocation(targets, euc2s);
    
    double3 euc1 = euc1s[0];
    double3 euc2 = euc2s[0];

    vector3 displacement(point3(2*euc1.x,2*euc1.y,2*euc1.z),point3(2*euc2.x,2*euc2.y,2*euc2.z));
    std::cout << euc1.x << ", " << euc1.y << ", " << euc1.z << "target:" << euc2.x << ", " <<euc2.y << ", " <<euc2.z << ", displacement:" << displacement << std::endl; 
    
     
    std::cout << source.faceIndex << ", " << source.x << std::endl;
    meshSpace->displaceParticle(source, displacement);
    std::cout << source.faceIndex << ", " << source.x << std::endl;
    }
