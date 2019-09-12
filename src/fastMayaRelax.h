#include <type_traits>
#include "fastRelax.h"
#include <maya/MMeshIntersector.h>
#include <maya/MFnMesh.h>
#include <maya/MObject.h>
#include <maya/MItGeometry.h>

// Double vs Float

typedef typename std::conditional<sizeof(FLOAT) == sizeof(double), MPoint, MFloatPoint>::type point_t;
typedef typename std::conditional<sizeof(FLOAT) == sizeof(double), MPointArray, MFloatPointArray>::type pointArray_t;

class MayaRelaxer: public Relaxer {
public:
	MMeshIntersector octree;
	MObject smoothMeshPar, smoothMesh;
	MayaRelaxer(MObject &mesh, MFnMesh &meshFn, MItGeometry& vertIter);

	void MayaRelaxer::quickRelax(
		MObject &mesh,
		const bool slide,
		const bool doReproject,
		const float taubinBias,
		const FLOAT iterations,
		FLOAT(*verts)[4]
	) const ;
	void buildOctree(MObject &mesh, bool slide, UINT divisions);
	void reprojectVerts(FLOAT(*verts)[4]) const;
	void reorderVerts(MObject &mesh, MFnMesh &meshFn, FLOAT(*reoVerts)[4]) const;
}

