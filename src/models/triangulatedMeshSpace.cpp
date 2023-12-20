#include "triangulatedMeshSpace.h"
/*! \file triangulatedMeshSpace.cpp */
#include <stdexcept>


void triangulatedMeshSpace::loadMeshFromFile(std::string filename, bool verbose)
    {
    if(verbose)
        printf("loading from file %s\n",filename.c_str());
    if(!CGAL::IO::read_polygon_mesh(filename, surface) || !CGAL::is_triangle_mesh(surface))
        {
        std::cerr << "Invalid input file." << std::endl;
        throw std::exception();
        };
    };

void triangulatedMeshSpace::geodesicDistance(meshPosition &p1, std::vector<meshPosition> &p2, std::vector<double> distances, std::vector<vector3> &startPathTangent, std::vector<vector3> &endPathTangent)
    {
    };

void triangulatedMeshSpace::displaceParticle(meshPosition &pos, vector3 &displacementVector)
    {
    pos.x += displacementVector;
    };
