#include "std_include.h"
#include <tclap/CmdLine.h>


#include "profiler.h"
#include "noiseSource.h"
#include "triangulatedMeshSpace.h"
#include "simulation.h"
#include "gradientDescent.h"
#include "velocityVerletNVE.h"
#include "simpleModel.h"
#include "gaussianRepulsion.h"
#include "harmonicRepulsion.h"
#include "vectorValueDatabase.h"
#include "cellListNeighborStructure.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


polyscope::SurfaceMesh *psMesh;
polyscope::PointCloud *psCloud;
std::vector<glm::vec3> pointCloud;
polyscope::CurveNetwork *psCurve;

char buf[256];
float sliderParameter1 = 1.;
int particleNumber = 50;
float deltaT = 0.01;
float interactionRange;
int timesteps = 100;
int stepsPerFrame = 10;
int index1 = 1;
int index2 = 2;
std::string loadMeshName;
bool configPopup = false;
bool dynPopup = false;
bool visPopup = false;
std::string basePathName = "geodesicPath";

std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

shared_ptr<simpleModel> configuration;
shared_ptr<cellListNeighborStructure> cellList;
shared_ptr<simulation> simulator;
shared_ptr<velocityVerletNVE> nve;
shared_ptr<gradientDescent> energyMinimizer;
shared_ptr<harmonicRepulsion> pairwiseForce; 
shared_ptr<triangulatedMeshSpace> meshSpace;

void getPointCloud(shared_ptr<simpleModel> model, vector<glm::vec3> &pos)
    {
    int N = model->N;
    model->fillEuclideanLocations();
    if(pointCloud.size()!= N)
        pointCloud.resize(configuration->N);
    for (int ii = 0; ii < N; ++ii)
        {
        double3 p = model->euclideanLocations[ii];
        pos[ii] = glm::vec3(p.x,p.y,p.z);
        }
    };
void setInteractionRange()
    {
    pairwiseForce = make_shared<harmonicRepulsion>(1.0,interactionRange);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);

    cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,interactionRange);
    configuration->setNeighborStructure(cellList);
    simulator->clearForceComputers();
    simulator->addForce(pairwiseForce);

    }

void setParticlesAndVelocities()
    {
    noiseSource noise(false);
    configuration=make_shared<simpleModel>(particleNumber);
    configuration->setSpace(meshSpace);
    configuration->setRandomParticlePositions(noise);
    configuration->setMaxwellBoltzmannVelocities(noise,sliderParameter1);
    cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,interactionRange);
    configuration->setNeighborStructure(cellList);
    pairwiseForce->setModel(configuration);

    simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);

    nve->setModel(configuration);
    simulator->addUpdater(nve,configuration);

    getPointCloud(configuration,pointCloud);
    psCloud = polyscope::registerPointCloud("particle positions", pointCloud);
    psCloud->setPointRadius(0.008);
    psCloud->addVectorQuantity("particle velocities",configuration->velocities);
    };

void runTimesteps()
    {
    polyscope::unshow();
    for (int ii = 1; ii <= timesteps; ++ii)
        {
        simulator->performTimestep();
        if(ii % stepsPerFrame ==0)
            {
            getPointCloud(configuration,pointCloud);
            psCloud->updatePointPositions(pointCloud);
            psCloud->addVectorQuantity("particle velocities",configuration->velocities);
            polyscope::frameTick();
            }
        }
    polyscope::show();
    };

void setNVEUpdater()
    {
    simulator->updaters.clear();
    nve=make_shared<velocityVerletNVE>(deltaT);
    nve->setModel(configuration);
    simulator->addUpdater(nve,configuration);
    };

void setGradientDescent()
    {
    simulator->updaters.clear();
    energyMinimizer=make_shared<gradientDescent>(deltaT);
    energyMinimizer->setModel(configuration);
    simulator->addUpdater(energyMinimizer,configuration);
    };

void loadNewMesh()
    {
    loadMeshName=buf;
    meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(loadMeshName,true);
    meshSpace->useSubmeshingRoutines(true,interactionRange,true);

    mesh.reset();
    geometry.reset();
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(loadMeshName);

    psMesh = polyscope::registerSurfaceMesh(
                    "Surface",
                    geometry->inputVertexPositions, mesh->getFaceVertexList(),
                    polyscopePermutations(*mesh));
    polyscope::view::resetCameraToHomeView();
    psMesh->setSmoothShade(true);
    setParticlesAndVelocities();
    };

void drawGeodesic()
    {
    shared_ptr<surfaceMeshShortestPath> smsp = make_shared<surfaceMeshShortestPath>(meshSpace->surface);
    meshPosition p1 = configuration->positions[index1];
    smspFaceLocation sourcePoint = meshPositionToFaceLocation(p1);
    AABB_tree localTree;
    smsp->build_aabb_tree(localTree);
    smsp->add_source_point(sourcePoint.first,sourcePoint.second);
    smsp->build_sequence_tree();
    meshPosition p2 = configuration->positions[index2];
    smspFaceLocation targetPoint = meshPositionToFaceLocation(p2);
    std::vector<point3> pathPoints;
    shortestPathResult geodesic = smsp->shortest_path_points_to_source_points(targetPoint.first, targetPoint.second,  std::back_inserter(pathPoints));
    psCurve = polyscope::registerCurveNetworkLine("geodesicPath", pathPoints);
    };

void drawGeodesics()
    {
    shared_ptr<surfaceMeshShortestPath> smsp = make_shared<surfaceMeshShortestPath>(meshSpace->surface);
    int N = configuration->getNumberOfParticles();
    for (int ii= 0; ii < N; ++ii)
        {
        meshPosition p1 = configuration->positions[ii];
        smspFaceLocation sourcePoint = meshPositionToFaceLocation(p1);
        AABB_tree localTree;
        smsp->build_aabb_tree(localTree);
        smsp->add_source_point(sourcePoint.first,sourcePoint.second);
        smsp->build_sequence_tree();
        int numberOfNeighbors = configuration->neighbors[ii].size();
        for (int jj =0; jj < numberOfNeighbors; ++jj)
            {
            int idx2 = configuration->neighbors[ii][jj];
            meshPosition p2 = configuration->positions[idx2];
            smspFaceLocation targetPoint = meshPositionToFaceLocation(p2);
            std::vector<point3> pathPoints;
            shortestPathResult geodesic = smsp->shortest_path_points_to_source_points(targetPoint.first, targetPoint.second,  std::back_inserter(pathPoints));
            string cnlName=basePathName+std::to_string(ii);
            psCurve = polyscope::registerCurveNetworkLine(cnlName, pathPoints);
            psCurve->setRadius(0.002);
            glm::vec3 newColor(0.867,0.289,0.285);
            psCurve->setColor(newColor);
            }
        }
    };

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback()
    {
    if (ImGui::Button("set configuration"))
        configPopup = true;
    if(configPopup)
        {
        if (ImGui::Begin("set configuration"))
            {
            ImGui::InputText("mesh filename", buf, sizeof(buf));         
            if (ImGui::Button("load new mesh"))
                {
                loadNewMesh();
                }
            ImGui::PushItemWidth(100);
            ImGui::InputInt("num points", &particleNumber);         
            ImGui::PopItemWidth();
            ImGui::SameLine();
            ImGui::PushItemWidth(100);
            ImGui::InputFloat("Temperature", &sliderParameter1);
            ImGui::PopItemWidth();
            if (ImGui::Button("set particle positions and velocities"))
                {
                setParticlesAndVelocities();
                }

            if(ImGui::Button("finished configuring",ImVec2(200, 0)))
                {
                configPopup= false;
                }
            ImGui::End();
            }
        }
    ImGui::SameLine();
    if (ImGui::Button("set dynamics"))
        dynPopup = true;
    if(dynPopup)
        {
        if (ImGui::Begin("choose dynamics"))
            {
            ImGui::InputFloat("delta t", &deltaT);
            ImGui::InputFloat("sigma", &interactionRange);
            if (ImGui::Button("set potential range"))
                {
                setInteractionRange();
                }
            if (ImGui::Button("set nve particle dynamics"))
                {
                setParticlesAndVelocities();
                setNVEUpdater();
                }
            if (ImGui::Button("set gradient descent"))
                {
                setGradientDescent();
                }

            if(ImGui::Button("finished setting dynamics",ImVec2(200, 0)))
                {
                dynPopup= false;
                }
            ImGui::End();
            }
        }
    ImGui::SameLine();
    if (ImGui::Button("visualize paths"))
        visPopup = true;
    if(visPopup)
        {
        if (ImGui::Begin("visualization"))
            {
            // ImGui::PushItemWidth(100);
            // ImGui::InputInt("particle a", &index1);
            // ImGui::SameLine();
            // ImGui::InputInt("particle b", &index2);
            // ImGui::PopItemWidth();
            if (ImGui::Button("draw interactions"))
                {
                drawGeodesics();
                }
            ImGui::SameLine();
            if (ImGui::Button("hide interactions"))
                {
                for (int ii = 0; ii < configuration->getNumberOfParticles(); ++ii)
                    polyscope::removeStructure(basePathName + std::to_string(ii));
                }

            if(ImGui::Button("finished visualization of paths",ImVec2(200, 0)))
            {
                visPopup= false;
                }
            ImGui::End();
            }
        }
    ImGui::Spacing();    

    ImGui::PushItemWidth(150);
    ImGui::InputInt("time steps", &timesteps);         
    ImGui::SameLine();
    ImGui::InputInt("frameUpdate", &stepsPerFrame);         
    ImGui::PopItemWidth();
    if (ImGui::Button("perform timesteps"))
        {
        runTimesteps();
        }

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
    ValueArg<int> iterationsArg("i","iterations","number of performTimestep calls to make",false,100,"int",cmd);
    ValueArg<int> saveFrequencyArg("s","saveFrequency","how often a file gets updated",false,100,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    ValueArg<double> interactionRangeArg("a","interactionRange","range ofthe interaction to set for both potential and cell list",false,1.,"double",cmd);
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
    particleNumber = N;
    int maximumIterations = iterationsArg.getValue();
    int saveFrequency = saveFrequencyArg.getValue();
    string meshName = meshSwitchArg.getValue();
    loadMeshName=meshName;
    std::strncpy(buf,loadMeshName.c_str(), sizeof(buf)-1);
    double dt = deltaTArg.getValue();
    double maximumInteractionRange= interactionRangeArg.getValue();
    interactionRange = maximumInteractionRange;
    double temperature = temperatureArg.getValue();
    bool verbose= verboseSwitch.getValue();
    bool reproducible = reproducibleSwitch.getValue();
    bool dangerous = dangerousSwitch.getValue(); //not used right now
    timesteps = maximumIterations;
    deltaT = dt;

    meshSpace=make_shared<triangulatedMeshSpace>();
    meshSpace->loadMeshFromFile(meshName,verbose);
    meshSpace->useSubmeshingRoutines(false);
    if(programBranch >=0)
        meshSpace->useSubmeshingRoutines(true,maximumInteractionRange,dangerous);

    configuration=make_shared<simpleModel>(N);
    configuration->setVerbose(verbose);
    configuration->setSpace(meshSpace);

    //set up the cellListNeighborStructure, which needs to know how large the mesh is
    cellList = make_shared<cellListNeighborStructure>(meshSpace->minVertexPosition,meshSpace->maxVertexPosition,maximumInteractionRange);
    configuration->setNeighborStructure(cellList);

    //for testing, just initialize particles randomly in a small space. Similarly, set random velocities in the tangent plane
    noiseSource noise(reproducible);
    configuration->setRandomParticlePositions(noise);
    configuration->setMaxwellBoltzmannVelocities(noise,temperature);

    //shared_ptr<gaussianRepulsion> pairwiseForce = make_shared<gaussianRepulsion>(1.0,.5);
    pairwiseForce = make_shared<harmonicRepulsion>(1.0,maximumInteractionRange);//stiffness and sigma. this is a monodisperse setting
    pairwiseForce->setModel(configuration);

    simulator=make_shared<simulation>();
    simulator->setConfiguration(configuration);
    simulator->addForce(pairwiseForce);

    nve=make_shared<velocityVerletNVE>(deltaT);
    nve->setModel(configuration);
    simulator->addUpdater(nve,configuration);

    for (int ii = 0; ii < 1; ++ii)
        simulator->performTimestep();

    
    // Initialize polyscope
    polyscope::init();
    // Set the callback function
    polyscope::state::userCallback = myCallback;
    // Load mesh
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshName);

    psMesh = polyscope::registerSurfaceMesh(
                    "Surface",
                    geometry->inputVertexPositions, mesh->getFaceVertexList(),
                    polyscopePermutations(*mesh));

    // Give control to the polyscope gui
    polyscope::options::programName = "curvedSpaceSimulations";
    polyscope::view::resetCameraToHomeView();
    psMesh->setSmoothShade(true);

    getPointCloud(configuration,pointCloud);
    psCloud = polyscope::registerPointCloud("particle positions", pointCloud);
    psCloud->addVectorQuantity("particle velocities",configuration->velocities);
    polyscope::show();
    return 0;
    };
