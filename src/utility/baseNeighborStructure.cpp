#include "baseNeighborStructure.h"
#include "std_include.h"

void baseNeighborStructure::initialize(std::vector<meshPosition> &_particles)
    {
    particles = _particles;
    indices.resize(particles.size());
    for(int ii = 0; ii < particles.size(); ++ii)
        {
        indices[ii]=ii;
        }
    }
/*!
Uses the "offset" to correctly identify "all particles other than some given index". For standard models there is no offset, but for MPI models the offset is important. The return value of this function is VERYLARGEDOUBLE, which should never be used for anything.
 */
double baseNeighborStructure::constructCandidateNeighborList(meshPosition &p, int particleIndex, std::vector<int> &candidateNeighborIndices, std::vector<meshPosition> &candidateParticles, int offset)
    {
    int N = particles.size();
    candidateNeighborIndices.resize(N-1);
    candidateParticles.resize(N-1);
    int elements = 0;
    for (int ii = 0; ii < N; ++ii)
        {
        if( (particleIndex + offset)!= ii)
            {
            candidateParticles[elements] = particles[ii];
            candidateNeighborIndices[elements] = indices[ii];
            elements +=1;
            };
        };
    return VERYLARGEDOUBLE;
    };
