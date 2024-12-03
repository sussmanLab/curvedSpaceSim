#ifndef openMeshSpace_H
#define openMeshSpace_H

#include "triangulatedMeshSpace.h"

class openMeshSpace : public triangulatedMeshSpace
    {
    public: 
        openMeshSpace(){};

	virtual void transportParticleAndVectors(meshPosition &pos, vector3 &displacementVector, vector<vector3> &transportVectors);

    protected:
        //!Implement boundary conditions; these functions will force particles to stop if they hit a boundary of the space
        virtual void boundaryEdge()
	    {
	    ERRORERROR("Boundary edge: base open mesh space does not have boundary conditions.");
	    }
        virtual void boundaryVertex()
	    {
	    ERRORERROR("Boundary vertex; base open mesh space does not have boundary conditions.");
	    };

	//!helper functions for boundary vertex behaviors -- these are designed to extend to more complex boundary conditions
        virtual void projectVectorsForBoundaryVertex(vector3 heading, vertexIndex interesctedV, vertexIndex boundaryV, vector3 fNormal, faceIndex face, vector<vector3> &transportVectors);
        virtual void getBoundaryVertexHeading(vertexIndex v, faceIndex &sourceFace, vector3 &displacementVector, pmpBarycentricCoordinates &sourceBarycentricLocation, vector<point3> vertexPositions, vector3 sourceNormal, std::pair<vector3, vertexIndex> &result);

    };

#endif
