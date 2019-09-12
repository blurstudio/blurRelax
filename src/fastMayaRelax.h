#include <type_traits>
#include <maya/MMeshIntersector.h>
#include <maya/MFnMesh.h>
#include <maya/MObject.h>
#include <maya/MItGeometry.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include "fastRelax.h"

// Double vs Float

typedef typename std::conditional<sizeof(FLOAT) == sizeof(double), MPoint, MFloatPoint>::type point_t;
typedef typename std::conditional<sizeof(FLOAT) == sizeof(double), MPointArray, MFloatPointArray>::type pointArray_t;

class MayaRelaxer: public Relaxer {
public:
	MMeshIntersector octree;
	MObject smoothMeshPar, smoothMesh;

	static MayaRelaxer Create(MObject &mesh, MFnMesh &meshFn, MItGeometry& vertIter,
		short borderBehavior, short hardEdgeBehavior, short groupEdgeBehavior
	);

	MayaRelaxer(
		short borderBehavior,
		short hardEdgeBehavior,
		short groupEdgeBehavior,
		std::vector<std::vector<size_t>> rawNeighbors,
		std::vector<std::vector<UCHAR>> rawHardEdges,
		const std::vector<UCHAR> &rawVertData
	) : Relaxer(borderBehavior, hardEdgeBehavior, groupEdgeBehavior, rawNeighbors, rawHardEdges, rawVertData)
	{ }
		
	void MayaRelaxer::quickRelax(
		MObject &mesh,
		const bool slide,
		const bool doReproject,
		const float taubinBias,
		const FLOAT iterations,
		FLOAT(*verts)[4]
	);
	void buildOctree(MObject &mesh, bool slide, UINT divisions);
	void reprojectVerts(FLOAT(*verts)[4]) const;
	void reorderVerts(MObject &mesh, MFnMesh &meshFn, FLOAT(*reoVerts)[4]) const;
};

