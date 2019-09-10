/*
MIT License

Copyright(c) 2018 Blur Studio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <maya/MTypeId.h> 
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MArrayDataHandle.h>
#include <maya/MGlobal.h>
#include <maya/MTypes.h>

#include <maya/MMeshIntersector.h>
#include <maya/MPointArray.h>
#include <maya/MFloatPointArray.h>

#include <maya/MItGeometry.h>
#include <maya/MPxDeformerNode.h> 

#include <maya/MFnPlugin.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnNumericAttribute.h>

#include <maya/MPoint.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MVector.h>
#include <maya/MMatrix.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MObject.h>
#include <maya/MMeshSmoothOptions.h>
#include <maya/MItMeshedge.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>

#include "external/xxhash.h"

#define CHECKSTAT(m) if (!status) {status.perror(m); return status;};

#define BB_NONE   0
#define BB_PIN    1
#define BB_SLIDE  2

#define HB_NONE   0
#define HB_PIN    1
#define HB_SLIDE  2

#define GB_NONE   0
#define GB_PIN    1
#define GB_SLIDE  2

#define SA_LAPLACIAN 0
#define SA_TAUBIN 1

#define DEFORMER_NAME "BlurRelax"

// Double vs Float
#define float_t double
#define point_t MPoint
#define pointArray_t MPointArray


void edgeProject(
	const float_t basePoints[][4],
	const std::vector<size_t> &group,
	const std::vector<size_t> &invOrder,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<std::vector<bool>> &hardEdges,
	const std::vector<UINT> &creaseCount,
	float_t smoothPoints[][4]
);

void quickLaplacianSmooth(
	float_t verts2d[][4],
	const size_t numVerts,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<float_t> &valence,
	const std::vector<float_t> &shiftVal,
	const std::vector<bool> &pinPoints,
	const float_t taubinBias=1.0
);

class BlurRelax : public MPxDeformerNode {
	public:
		BlurRelax();
		~BlurRelax() override;

		static void* creator();
		static MStatus initialize();

		MStatus deform(MDataBlock& dataBlock, MItGeometry& vertIter, const MMatrix& matrix, UINT multiIndex) override;

	public:
		// local node attributes
		static MObject aIterations;
		static MObject aBorderBehavior;
		static MObject aHardEdgeBehavior;
		static MObject aGroupEdgeBehavior;
		static MObject aReproject;
		static MObject aTaubinBias;
		static MTypeId id;
	private:
		// Hash checking variables
		unsigned long long faceHash = 0;
		unsigned long long edgeHash = 0;
		unsigned long long groupHash = 0;
		unsigned long long smoothHash = 0;
		short bbCheck = 255;
		short hbCheck = 255;
		short gbCheck = 255;

		// storage for this data doesn't change unless the hashes do
		std::vector<size_t> order;
		std::vector<size_t> invOrder;
		std::vector<std::vector<size_t>> neighbors;
		std::vector<std::vector<bool>> hardEdges;
		std::vector<float_t> shiftVal; // normally 0.5; but it's 0.25 if on a hard edge
		std::vector<float_t> valence; // as float for vectorizing
		std::vector<bool> pinPoints;
		std::vector<UINT> creaseCount;


		MStatus getTrueWeights(
				MObject &mesh,
				MDataBlock &dataBlock,
				UINT index,
				std::vector<float> &weightVals,
				float envelope
				) const;

		void BlurRelax::buildQuickData(
			MObject &mesh,
			MItGeometry& vertIter,
			short borderBehavior,
			short hardEdgeBehavior,
			short groupEdgeBehavior,
			bool reproject,
			std::vector<char> &group
		);

		void quickRelax(
			MObject &mesh,
			const short borderBehavior,
			const short hardEdgeBehavior,
			const short groupEdgeBehavior,
			const bool reproject,
			const float taubinBias,
			const float_t iterations,
			const UINT numVerts,
			const std::vector<char> &group,
			float_t(*verts)[4] // already resized
		);
};



