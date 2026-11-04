#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MItGeometry.h>
#include <maya/MMatrixArray.h>
#include <maya/MMeshIntersector.h>
#include <maya/MObject.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>

#include <type_traits>

#include "fastRelax.h"

// Double vs Float

typedef
    typename std::conditional<sizeof(FLOAT) == sizeof(double), MPoint, MFloatPoint>::type point_t;
typedef
    typename std::conditional<sizeof(FLOAT) == sizeof(double), MPointArray, MFloatPointArray>::type
        pointArray_t;

class MayaRelaxer : public Relaxer {
  public:
    static MayaRelaxer *Create(
        MObject &mesh, MFnMesh &meshFn, MItGeometry &vertIter, short borderBehavior,
        short hardEdgeBehavior, short groupEdgeBehavior
    );

    MayaRelaxer(
        short borderBehavior, short hardEdgeBehavior, short groupEdgeBehavior,
        std::vector<std::vector<size_t>> rawNeighbors, std::vector<std::vector<UCHAR>> rawHardEdges,
        const std::vector<UCHAR> &rawVertData
    )
        : Relaxer(
              borderBehavior, hardEdgeBehavior, groupEdgeBehavior, rawNeighbors, rawHardEdges,
              rawVertData
          )
    {
    }

    pointArray_t quickRelax(
        MObject &mesh, const bool slide, const bool doReproject, const UINT reprojectDivisons,
        const float taubinBias, const FLOAT iterations
    );

    // Relax an arbitrary set of per-vertex values using this relaxer's topology.
    // The mesh is only used if doReproject is true
    pointArray_t quickRelax(
        const pointArray_t &points, MObject &mesh, const bool slide, const bool doReproject,
        const UINT reprojectDivisons, const float taubinBias, const FLOAT iterations
    );

    void buildOctree(MObject &mesh, bool slide, UINT divisions);
    void reprojectVerts(FLOAT (*verts)[NUM_COMPS], MMeshIntersector &octree) const;

    void reorderVerts(MObject &mesh, MFnMesh &meshFn, FLOAT (*reoVerts)[NUM_COMPS]) const;
    void reorderVerts(const pointArray_t &mpa, FLOAT (*reoVerts)[NUM_COMPS]) const;

    void revertVerts(FLOAT (*reoVerts)[NUM_COMPS], FLOAT (*reverted)[NUM_COMPS]) const;
    pointArray_t revertVerts(FLOAT (*reoVerts)[NUM_COMPS]) const;

    std::vector<MMatrix> getVertMatrices(MObject &mesh, MFnMesh &meshFn) const;
};
