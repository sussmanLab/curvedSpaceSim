#include "simpleModel.h"
/*! \file simpleModel.cpp */

/*!
 * Set the size of basic data structures...
*/
simpleModel::simpleModel(int n) :
    N(n)
    {
    cout << "initializing a model with "<< N << " particles" << endl;
    initializeSimpleModel(n);
    };

/*!
 * actually set the array sizes. positions, velocities, forces are zero
 * masses are set to unity
*/
void simpleModel::initializeSimpleModel(int n)
    {
    N=n;

    positions.resize(n);
    euclideanLocations.resize(n);
    velocities.resize(n,vector3(0.,0.,0.));
    forces.resize(n,vector3(0.,0.,0.));
    masses.resize(n);
    masses = vector<double>(n,1.0);
    types.resize(n,0);

    neighbors.resize(n);
    neighborVectors.resize(n);
    neighborDistances.resize(n);
    //radii.resize(n);

    neighborStructure = make_shared<baseNeighborStructure>();
    };

void simpleModel::fillEuclideanLocations()
    {
    space->meshPositionToEuclideanLocation(positions,euclideanLocations);
    }

void simpleModel::moveParticles(vector<vector3> &disp)
    {
    if(!particleShiftsRequireVelocityTransport)
        {
        for(int ii = 0; ii < N; ++ii)
	    {
	    //cout << "displacing particle " << ii << endl; //debug statement to make sure shifts are finishing
            vector3 fModify = forces[ii];
            space->displaceParticle(positions[ii],disp[ii], fModify);
	    forces[ii] = fModify;
	    }
        }
    else 
	{
        for(int ii = 0; ii < N; ++ii)
	    {
	    vector3 fModify = forces[ii];
	    space->transportParticleAndVelocity(positions[ii],velocities[ii],disp[ii], forces[ii]);
            forces[ii] = fModify;
	    }
	}
    };

void simpleModel::findNeighbors(double maximumInteractionRange)
    {
    //first, determine any needed neighbor structure initialization
    neighborStructure->setInteractionRange(maximumInteractionRange);
    bool euclideanNeighborsMeshPositions = (neighborStructure->requireEuclideanPositions && !space->positionsAreEuclidean);
    if(euclideanNeighborsMeshPositions)
        {
        space->meshPositionToEuclideanLocation(positions,euclideanMeshPosition);
        neighborStructure->initialize(euclideanMeshPosition);
        }
    else
        neighborStructure->initialize(positions);

    double meanNN = 0;
    for (int ii =0; ii < N; ++ii)
        {
        //use the neighborStructure to find candidate neighbors. This function fills the indices of the neighbors[ii] data structure and populate the positions of the corresponding targetParticles
        double largestNeighborDistance;
        vector<meshPosition> targetParticles;
        //if needed, target a euclidean position
        if(euclideanNeighborsMeshPositions)
            largestNeighborDistance = neighborStructure->constructCandidateNeighborList(euclideanMeshPosition[ii], ii, neighbors[ii], targetParticles);
        else
            largestNeighborDistance = neighborStructure->constructCandidateNeighborList(positions[ii], ii, neighbors[ii], targetParticles);

        //if needed, re-fill targetParticles with space-appropriate data
        if(euclideanNeighborsMeshPositions)
            {
            for(int jj = 0; jj < neighbors[ii].size();++jj)
                targetParticles[jj] = positions[neighbors[ii][jj]];
            };
        //additionally use the space to populate the list of distances and separation vectors
        vector<vector3> placeholderVector;
        vector<vector3> tangentVector;
        vector<double> distances;
        space->distance(positions[ii],targetParticles,distances,tangentVector,placeholderVector,largestNeighborDistance);
        neighborDistances[ii] = distances;
        neighborVectors[ii] = tangentVector;
        meanNN += distances.size();
        };
    if(verbose)
        printf("mean number of neighbors = %f\n",meanNN/N);
    };

void simpleModel::findBoundaryParticles(vector<int> &result)
    {
    result.reserve(positions.size()/2); //way more than it needs
    for (int ii = 0; ii < positions.size(); ii++)
        {
        faceIndex iiFace = faceIndex(positions[ii].faceIndex);
        //is_border only works with halfedges, so we need to dig out the associated edges
        halfedgeIndex hf(tMeshSpace->surface.halfedge(iiFace));
        for(halfedgeIndex hi : halfedges_around_face(hf, tMeshSpace->surface))
            {
            edgeIndex ei = tMeshSpace->surface.edge(hi);
            if(tMeshSpace->surface.is_border(ei))
                {
                result.push_back(ii);
                break;
                }
            }
        }
    }

void simpleModel::clampBarycentricCoordinatesToFace(point3 &barycentricWeights)
    {
    double w1 = barycentricWeights.x();
    double w2 = barycentricWeights.y();
    double w3 = barycentricWeights.z();
    double tol = clampTolerance;
    if ((w1 < -tol) || (w2 < -tol) || (w3 < -tol)) cout << "While clamping, weight negative in its face, weights " << w1 << ", " << w2 << ", " << w3 << endl;
    if (abs(w1) < tol) w1 = tol;
    if (abs(w2) < tol) w2 = tol;
    if (abs(w3) < tol) w3 = tol;

    w1 = w1/(w1+w2+w3);
    w2 = w2/(w1+w2+w3);
    w3 = w3/(w1+w2+w3);

    barycentricWeights = point3(w1,w2,w3);
    }

void simpleModel::R3PositionsToMeshPositions(triangleMesh &mesh, vector<point3> r3positions, vector<meshPosition> &simPositions)
    {
    // This function converts a set of R3 positions, expressed as a vector of Point3 objects, to mesh positions.
    // BE CAUTIOUS -- because R3 positions are generally very slightly off-mesh, we use CGAL's locate functions
    // to find the nearest position on the mesh. The function will *not* throw an error if you use points that
    // are very far away from the mesh, even though the resulting mesh positions will likely be drastically different
    // from the points you fed in.
    AABB_tree tree;
    PMP::build_AABB_tree(mesh, tree);

    for (point3 pos : r3positions)
        {
        pmpFaceLocation locateOutput = PMP::locate_with_AABB_tree(pos, tree, mesh);
        point3 baryWeights = point3(locateOutput.second[0],locateOutput.second[1],locateOutput.second[2]);
        clampBarycentricCoordinatesToFace(baryWeights);
        int face = locateOutput.first;
        simPositions.push_back(meshPosition(baryWeights, face));
        }
    }

void simpleModel::setMeshPositionsFromR3File(string filename,  triangleMesh &mesh)
    {
    vector<meshPosition> simPositions;
    //simPositions needs to be set to the same length as the number of positions in input file
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Failed to open position file." << std::endl;
    }

    vector<point3> points;

    // Read the file line by line
    while (getline(file, line))
        {
        stringstream ss(line);
        string entry;
        vector<double> entries;
        point3 point;

        // Split the line by comma and read each entry
        while (getline(ss, entry, ','))
            {
            double value;
            istringstream(entry) >> value; //allows the string to convert to double
            entries.push_back(value);
            }

        // Check if the line contains enough entries
        if (entries.size() != 3)
            {
            cerr << "Error: input file had invalid number of entries on a line. Skipping line." << endl;
            continue; // Skip this line and move to the next one
            }

        // Use entries to make point3
        point = point3(entries[0], entries[1], entries[2]);

        // Add the tuple to the vector
        points.push_back(point);
        }
    simPositions.reserve(points.size());
    R3PositionsToMeshPositions(mesh, points, simPositions);
    setParticlePositions(simPositions);
    }


void simpleModel::setParticlePositions(vector<meshPosition> &newPositions)
    {
    if(N !=newPositions.size())
        initializeSimpleModel(newPositions.size());
    for (int pp = 0;pp < N; ++pp)
        {
        positions[pp] = newPositions[pp];
        };
    };

void simpleModel::setRandomParticlePositions(noiseSource &noise)
    {
    for(int pp = 0; pp < N; ++pp)
        space->randomPosition(positions[pp], noise);
    }

/*!
Assumes all particles have unit mass (for now)
*/
void simpleModel::setMaxwellBoltzmannVelocities(noiseSource &noise, double T)
    {
    for (int pp = 0; pp < N; ++pp)
        {
        space->randomVectorAtPosition(positions[pp], velocities[pp],noise);
        velocities[pp] *= sqrt(T);
        }
    };

/*!
It is possible that in some complex models it may be easier to implement a force
calculation from within the model (see, e.g., the way cellGPU handles
cell-model-like forces).  In that case, the simulation should call the model's
internal force function.  This is currently not implement in the simulation
classes, but would be easy to do.
TODO
*/
void simpleModel::computeForces(bool zeroOutForces)
    {
    if(zeroOutForces)
        {
        for(int ii = 0; ii < N; ++ii)
            {
            forces[ii] = vector3(0.0,0.0,0.0);
            };
        }
    };

