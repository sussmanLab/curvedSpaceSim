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

    noiseSource noise(reproducible);



    shared_ptr<triangulatedMeshSpace> meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);

    int nFaces = meshSpace->surface.number_of_faces();
    printf("%i\n",meshSpace->surface.number_of_vertices());

    int randomFace = noise.getInt(0,nFaces-1);
    faceDescriptor fdTest(randomFace);
    double3 bary= noise.getRandomBarycentricSet();
    pmpBarycentricCcoordinates baryTest = {bary.x,bary.y,bary.z};
    pmpFaceLocation randomLocationOnRandomFace(randomFace,baryTest);
    
    
    smspFaceLocation smspRL = randomLocationOnRandomFace;


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

    return 0;
    };
