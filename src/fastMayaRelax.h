#include <type_traits>
#include <maya/MMeshIntersector.h>
#include <maya/MFnMesh.h>
#include <maya/MObject.h>
#include <maya/MItGeometry.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MMatrixArray.h>
#include "fastRelax.h"

// Double vs Float

typedef typename std::conditional<sizeof(FLOAT) == sizeof(double), MPoint, MFloatPoint>::type point_t;
typedef typename std::conditional<sizeof(FLOAT) == sizeof(double), MPointArray, MFloatPointArray>::type pointArray_t;

class MayaRelaxer: public Relaxer {
public:
	std::vector<int> connFaces;

	static MayaRelaxer* Create(MObject &mesh, MFnMesh &meshFn, MItGeometry& vertIter,
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
		
	pointArray_t quickRelax(
		MObject &mesh,
		const bool slide,
		const bool doReproject,
		const UINT reprojectDivisons,
		const float taubinBias,
		const FLOAT iterations
	);

	void buildOctree(MObject &mesh, bool slide, UINT divisions);
	void reprojectVerts(FLOAT(*verts)[NUM_COMPS], MMeshIntersector &octree) const;

	void reorderVerts(MObject &mesh, MFnMesh &meshFn, FLOAT(*reoVerts)[NUM_COMPS]) const;
	void reorderVerts(pointArray_t mpa, FLOAT(*reoVerts)[NUM_COMPS]) const;

	void revertVerts(FLOAT(*reoVerts)[NUM_COMPS], FLOAT(*reverted)[NUM_COMPS]) const;
	pointArray_t revertVerts(FLOAT(*reoVerts)[NUM_COMPS]) const;


	MMatrix getMatrixAtPoint(MObject &mesh, MFnMesh &meshFn, MItMeshVertex &vertIt) const;
	std::vector<MMatrix> getVertMatrices(MObject &mesh, MFnMesh &meshFn) const;

};

