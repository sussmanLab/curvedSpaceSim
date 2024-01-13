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

    noiseSource noise(reproducible);



    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);

    int nFaces = meshSpace->surface.number_of_faces();
    printf("%i\n",meshSpace->surface.number_of_vertices());

    int randomFace = noise.getInt(0,nFaces-1);
    faceDescriptor fdTest(randomFace);
    double3 bary= noise.getRandomBarycentricSet();
    pmpBarycentricCoordinates baryTest = {bary.x,bary.y,bary.z};
    pmpFaceLocation randomLocationOnRandomFace(randomFace,baryTest);


    smspFaceLocation smspRL = randomLocationOnRandomFace;
profiler p1("geodesic");
profiler p2("shift");

vector<meshPosition> pos2;
vector<pmpFaceLocation> faceTargets;
for(int ii = 0; ii < 50; ++ii)
    {
    smspBarycentricCoordinates  target;
    target[0]=noise.getRealUniform(-2.5,2.5);
    target[1]=noise.getRealUniform(-2.5,2.5);
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
double maxDist = .8;//just a random scale for testing the submesher
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
