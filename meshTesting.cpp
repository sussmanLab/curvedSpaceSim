#include "std_include.h"
#include <tclap/CmdLine.h>

#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "euclideanSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "cellListNeighborStructure.h"
#include "meshUtilities.h"
#include "submesher.h"

using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in curved space!",' ',"V0.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> particleNumberSwitchArg("n","number","number of particles to simulate",false,20,"int",cmd);
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,1000,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    SwitchArg reproducibleSwitch("r","reproducible","reproducible random number generation", cmd, true);
    SwitchArg dangerousSwitch("d","dangerousMeshes","meshes where submeshes are dangerous", cmd, false);
    SwitchArg verboseSwitch("v","verbose","output more things to screen ", cmd, false);
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);

    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    int N = particleNumberSwitchArg.getValue();
    int maximumIterations = iterationsArg.getValue();
    string meshName = meshSwitchArg.getValue();
    double dt = 0.01;
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = dangerousSwitch.getValue();
    double maxDist = interactionRangeArg.getValue();

    if(dangerous && programBranch >= 0)
        printf("%i %f %i\n",N,dt,maximumIterations);
    noiseSource noise(reproducible);


if(programBranch ==1)
{
    //meshSpace might use dangerous submeshing, meshSpace2 won't
    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(true,maxDist,dangerous);
    shared_ptr<triangulatedMeshSpace> meshSpace2=make_shared<triangulatedMeshSpace>();
    meshSpace2->loadMeshFromFile(meshName,verbose);
    meshSpace2->useSubmeshingRoutines(true,maxDist,false);

    shared_ptr<simpleModel> configuration1=make_shared<simpleModel>(N);
    shared_ptr<simpleModel> configuration2=make_shared<simpleModel>(N);
    configuration1->setVerbose(verbose);
    configuration1->setSpace(meshSpace);
    configuration2->setVerbose(verbose);
    configuration2->setSpace(meshSpace2);

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    shared_ptr<cellListNeighborStructure> cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maxDist);
    shared_ptr<cellListNeighborStructure> cellList2 = make_shared<cellListNeighborStructure>(meshSpace2->minVertexPosition,meshSpace2->maxVertexPosition,maxDist);
    configuration1->setNeighborStructure(cellList);
    configuration2->setNeighborStructure(cellList2);
    noiseSource noise1(reproducible);
    noiseSource noise2(reproducible);
    configuration1->setRandomParticlePositions(noise1);
    configuration2->setRandomParticlePositions(noise2);
    shared_ptr<harmonicRepulsion> pairwiseForce1 = make_shared<harmonicRepulsion>(1.0,maxDist);//stiffness and sigma. this is a monodisperse setting
    shared_ptr<harmonicRepulsion> pairwiseForce2 = make_shared<harmonicRepulsion>(1.0,maxDist);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce1->setModel(configuration1);
    pairwiseForce2->setModel(configuration2);

    pairwiseForce1->computeForces(configuration1->forces);
    pairwiseForce2->computeForces(configuration2->forces);

    for(int ii = 0; ii < N; ++ii)
        {
        int n1 = 0;
        int n2 = 0;
        for(int jj = 0; jj < configuration1->neighborDistances[ii].size(); ++jj)
            if(configuration1->neighborDistances[ii][jj] < maxDist) n1 +=1;
        for(int jj = 0; jj < configuration2->neighborDistances[ii].size(); ++jj)
            if(configuration2->neighborDistances[ii][jj] < maxDist) n2 +=1;
        if(n1!=n2)
            {
            printf("i %i\t n1 %i n2 %i\n",ii,n1,n2);
            }
        }
}

if(programBranch ==0)
{
        shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
        meshSpace->loadMeshFromFile(meshName,verbose);

        int nFaces = meshSpace->surface.number_of_faces();
        printf("%i\n",meshSpace->surface.number_of_vertices());
        int randomFace = noise.getInt(0,nFaces-1);
        faceDescriptor fdTest(randomFace);
        double3 bary= noise.getRandomBarycentricSet();
        pmpBarycentricCoordinates baryTest = {bary.x,bary.y,bary.z};
        pmpFaceLocation randomLocationOnRandomFace(randomFace,baryTest);


        profiler p1("geodesic");
        profiler p2("shift");

        vector<meshPosition> pos2;
        vector<pmpFaceLocation> faceTargets;
        for(int ii = 0; ii < N; ++ii)
        {
                smspBarycentricCoordinates  target;
                target[0]=noise.getRealUniform(-4,4);
                target[1]=noise.getRealUniform(-4,4);
                target[2]=1-target[0]-target[1];
                pmpFaceLocation targetFL(fdTest, target);

                point3 targetPoint = PMP::construct_point(targetFL,meshSpace->surface);
                point3 sourcePoint = PMP::construct_point(randomLocationOnRandomFace,meshSpace->surface);
                vector3 displacementVector(sourcePoint,targetPoint);
                meshPosition testPoint;
                testPoint.x = point3(randomLocationOnRandomFace.second[0],randomLocationOnRandomFace.second[1],randomLocationOnRandomFace.second[2]);
                testPoint.faceIndex = fdTest;
                p2.start();
                meshSpace->displaceParticle(testPoint, displacementVector);
                p2.end();
                pos2.push_back(testPoint);
                targetFL.first=(faceIndex) testPoint.faceIndex;

                targetFL.second[0] = testPoint.x[0];
                targetFL.second[1] = testPoint.x[1];
                targetFL.second[2] = testPoint.x[2];

                faceTargets.push_back(targetFL);
        }

        vector<double> distances;
        vector<vector3> startPath;
        vector<vector3> endPath;
        meshPosition rlrf;
        rlrf.x = point3(randomLocationOnRandomFace.second[0],randomLocationOnRandomFace.second[1],randomLocationOnRandomFace.second[2]);
        rlrf.faceIndex = randomLocationOnRandomFace.first;
        p1.start();
        meshSpace->distance(rlrf,pos2,distances,startPath,endPath);
        p1.end();
        p1.print();
        p2.print();
        cout <<endl;
        vector<faceLocation> faceTargetsForSubmesh;
        vector<meshPosition> mpTargetsForSubmesh;
        for(int ii = 0; ii < distances.size();++ii)
        {
                if(distances[ii] < maxDist)
                {
                        point3 pp = PMP::construct_point(faceTargets[ii],meshSpace->surface);
                        meshPosition mp;
                        mp.x = point3(faceTargets[ii].second[0],faceTargets[ii].second[1],faceTargets[ii].second[2]);
                        mp.faceIndex = faceTargets[ii].first;
                        printf("{%f,%f,%f},",pp[0],pp[1],pp[2]);
                        faceTargetsForSubmesh.push_back(faceTargets[ii]);
                        mpTargetsForSubmesh.push_back(mp);
                }
        }
        cout <<endl <<" assessing performance on a submesh with " << faceTargetsForSubmesh.size() << " targets" << endl;
        submesher submeshAssistant;
        profiler pSubmesh("submesh construction");
        profiler pSubG("submesh construction and geodesic finding");
        pSubmesh.start();
        std::unordered_map<faceIndex,int> fMap;
        std::unordered_map<vertexIndex,int> vMap;
        triangleMesh submesh = submeshAssistant.constructSubmeshFromSourceAndTargets(meshSpace->surface, randomLocationOnRandomFace, faceTargetsForSubmesh,maxDist,vMap,fMap);
        pSubmesh.end();
        pSubmesh.print();

        vector<double> distancesSubmesh;
        profiler p0("geodesic subset");
        p0.start();
        meshSpace->distance(rlrf,mpTargetsForSubmesh,distances,startPath,endPath);
        p0.end();
        p0.print();
        meshSpace->useSubmeshingRoutines(true,maxDist,dangerous);

        pSubG.start();
        meshSpace->distance(rlrf,mpTargetsForSubmesh,distancesSubmesh,startPath,endPath);
        pSubG.end();
        pSubG.print();
        double jj = 0;
        for (int ii = 0; ii < distances.size(); ++ii)
                jj += distances[ii] - distancesSubmesh[ii];
        printf("total difference in computed distances between full and submesh routines: %f\n", jj);
}
//test the throughVertex routine with a known case
    if(programBranch == -1)
    {
    meshName = "../exampleMeshes/torus_isotropic_remesh.off"; 
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
    };
/*
//spot test of edge intersection detection
for(int ii = 0; ii < 10; ++ii)
{
randomFace = noise.getInt(0,nFaces-1);
faceDescriptor fDtest2(randomFace);
std::vector<point3> vertices;
getVertexPositionsFromFace(meshSpace->surface, faceDescriptor(randomFace), vertices);
printf("{");
printPoint(vertices[0]);
printPoint(vertices[1]);
printPoint(vertices[2]);
double3 sourcePoint= noise.getRandomBarycentricSet();
smspBarycentricCoordinates source, target, i1,i2,i3;
source[0] = sourcePoint.x;
source[1] = sourcePoint.y;
source[2] = sourcePoint.z;
target[0]=noise.getRealUniform(-2,2);
target[1]=noise.getRealUniform(-2,2);
target[2]=1-target[0]-target[1];
bool ii1 = intersectionBarycentricLinesV1V2(source,target,i1);
bool ii2 = intersectionBarycentricLinesV2V3(source,target,i2);
bool ii3 = intersectionBarycentricLinesV3V1(source,target,i3);
printBary(source);
printBary(target);
if (ii1) printBary(i1);
if (ii2) printBary(i2);
if (ii3)printBary(i3);
printf("}\n");
}
 */
/*
   smspFaceLocation smspRL = randomLocationOnRandomFace;
//spot testing of simple timing
point3 targetPoint = pathFinder.point(smspRL.first,smspRL.second);
profiler p2("smsp initialize");
p2.start();
surfaceMeshShortestPath pathFinder(meshSpace->surface);
p2.end();

profiler p3("point");
p3.start();
point3 targetPoint = pathFinder.point(smspRL.first,smspRL.second);
p3.end();

profiler p1("locate");
p1.start();
AABB_tree tree;
pathFinder.build_aabb_tree(tree);
smspFaceLocation loc2 = pathFinder.locate<AABB_face_graph_traits>(targetPoint,tree);
p1.end();

p1.print();
p2.print();
p3.print();

printf("(%f %f %f)\n",targetPoint[0],targetPoint[1],targetPoint[2]);
printf("(%f %f %f) (%f %f %f) %i %i\n", 
smspRL.second[0],smspRL.second[1],smspRL.second[2],
loc2.second[0],loc2.second[1],loc2.second[2],
(int)loc2.first,randomFace
);
profiler p4("getVP");
p4.start();
std::vector<point3> testV;
getVertexPositionsFromFace(meshSpace->surface, loc2.first,testV);
p4.end();
printf("(%f %f %f) (%f %f %f) (%f %f %f)\n ",
testV[0][0],testV[0][1],testV[0][2],
testV[1][0],testV[1][1],testV[1][2],
testV[2][0],testV[2][1],testV[2][2]
);
p4.print();
 */
return 0;
};
