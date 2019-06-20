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
) {
	/*
	This code takes the edges that have been defined as hard
	and reprojects them onto the original neighboring hard edges

	Perhaps in the future, I can project onto non-branching
	strings of edges. But for now, it's just direct neighbors
	*/

	std::vector<size_t> neigh;
	for (size_t gidx = 0; gidx < group.size(); ++gidx) {
		// If we have "hard edges" we have already removed
		// any soft edges, so if the crease count for this vertex is zero
		// we can just completely skip it. And if it's >0, then we can
		// assume that all stored neighbors are hard.
		size_t idx = invOrder[group[gidx]];

		float_t *avg = smoothPoints[idx];
		const float_t *basePos = basePoints[idx];
		if (creaseCount[idx] != 2) { continue; }
		float_t keep[3], delta[3], edge[3];

		neigh.clear();
		for (size_t i = 0; i < neighbors.size(); ++i) {
			if (idx < neighbors[i].size()) {
				neigh.push_back(neighbors[i][idx]);
			}
			else break;
		}

		// TODO: Find the hard edge strings and reproject onto those
		float_t minLen = std::numeric_limits<float_t>::max();
		bool found = false;

		// avg - basePos
		for (size_t x = 0; x < 3; ++x)
			delta[x] = avg[x] - basePos[x];

		for (size_t i = 0; i < neigh.size(); ++i) {
			size_t n = neigh[i];

			// normalized(prevPoints[n] - basePos)
			for (size_t x = 0; x < 3; ++x)
				edge[x] = basePoints[n][x] - basePos[x];
			float_t elen = sqrt(edge[0] * edge[0] + edge[1] * edge[1] + edge[2] * edge[2]);
			for (size_t x = 0; x < 3; ++x)
				edge[x] /= elen;

			// dot(delta, edge)
			float_t dd = delta[0] * edge[0] + delta[1] * edge[1] + delta[2] * edge[2];

			float_t dn[3] = { 0.0f, 0.0f, 0.0f };
			if (dd > 0.0) {
				for (size_t x = 0; x < 3; ++x)
					dn[x] = edge[x] * dd;
			}
			float_t xx[3];

			// delta - dn
			for (size_t x = 0; x < 3; ++x)
				xx[x] = delta[x] - dn[x];

			// dot(xx, xx)
			float_t len2 = xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2];

			if (len2 < minLen) {
				minLen = len2;
				for (size_t x = 0; x < 3; ++x)
					keep[x] = dn[x];
				found = true;
			}
		}

		if (found) {
			for (size_t x = 0; x < 3; ++x)
				avg[x] = basePos[x] + keep[x];
		}

	}
}

void quickLaplacianSmooth(
	float_t verts2d[][4],
	const size_t numVerts,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<float_t> &valence,
	const std::vector<float_t> &shiftVal,
	const std::vector<bool> &pinPoints,
	const float_t taubinBias=1.0
) {
	/*
	All the crazy hoops I've jumped through are to make auto-vectorization work

	neighbors and neighborOffsets work together to be a kind of pseudo-transpose
	of a std::vector of neighbors per vert with that vector sorted so the verts
	with the most neighbors were at the top.
	neighborOffsets contains the index offsets of vertices with at least [index] of neighbors
	*/

	// First, get verts as a single pointer to the contiguous memory stored in (verts2d*)[4]
	float_t* verts = &(verts2d[0][0]);

	// number of nonzero valence
	size_t nzv = neighbors[0].size();

	// number of nonzero components
	size_t nzc = 4 * nzv;

	// The __restrict keyword tells the compiler that *outComp
	// is not pointed to by any other pointer in this scope
	// This allows for auto-vectorization
	float_t * __restrict outComp = new float_t[nzc];
	memset(outComp, 0, nzc*sizeof(float_t));

	for (size_t ncIdx = 0; ncIdx < neighbors.size(); ++ncIdx) {
		const auto &nCol = neighbors[ncIdx];
		size_t nColCount = nCol.size();
		for (size_t i = 0; i < nColCount; ++i) {
			size_t nci = 4 * nCol[i];
			for (size_t j = 0; j < 4; ++j) {
				outComp[4*i+j] = outComp[4*i+j] + verts[nci+j];
			}
		}
	}

	for (size_t i = 0; i < nzc; ++i) {
		outComp[i] = shiftVal[i] * taubinBias * ((outComp[i] / valence[i]) - verts[i]) + verts[i];
	}

	memcpy(verts, outComp, nzc*sizeof(float_t));
	delete outComp;
}


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

MTypeId BlurRelax::id(0x001226FD);

// local attributes
MObject BlurRelax::aIterations;
MObject BlurRelax::aBorderBehavior;
MObject BlurRelax::aHardEdgeBehavior;
MObject BlurRelax::aGroupEdgeBehavior;
MObject BlurRelax::aReproject;
MObject BlurRelax::aTaubinBias;


BlurRelax::BlurRelax() {}
BlurRelax::~BlurRelax() {}

void* BlurRelax::creator() {
	return new BlurRelax();
}

MStatus BlurRelax::getTrueWeights(MObject &mesh, MDataBlock &dataBlock, UINT index, std::vector<float> &weightVals, float envelope) const {
	// Maya might do some tricky stuff like not store the weight array at all for certain weight
	// values so we can't count on an array existing in the weightList. For the OpenCL Kernel
	// we want an array with one weight in it per vertex, we need to build it carefully here.

	// Two possibilities: we could have a sparse array in weightList[multiIndex] or there could be nothing in weightList[multiIndex].
	// if nothing is there then all the weights at 1.0f.

	MStatus status;
	MFnMesh meshFn(mesh);
	UINT numVerts = meshFn.numVertices();
	weightVals.reserve(numVerts);

	// Get a handle to the weight array we want.
	MArrayDataHandle weightList = dataBlock.outputArrayValue(MPxDeformerNode::weightList, &status);
	if (!status) return status; // we should always be able to get a weightList
	status = weightList.jumpToElement(index);
	// it is possible that the jumpToElement fails. In that case use the envelope value
	if (!status) {
		weightVals.resize(numVerts, envelope);
	}
	else {
		MDataHandle hWeightsStructure = weightList.inputValue(&status);
		if (!status) return status;
		MArrayDataHandle hWeights = hWeightsStructure.child(MPxDeformerNode::weights);
		if (!status) return status;

		// number of non-zero weights
		unsigned int numWeights = hWeights.elementCount(&status);
		if (!status) return status;

		// we're building a list with a weight per vertex, even if the weight is zero
		unsigned int weightIndex = 0;
		for (unsigned int i = 0; i < numWeights; i++, hWeights.next()) {
			unsigned int weightsElementIndex = hWeights.elementIndex(&status);
			while (weightIndex < weightsElementIndex) {
				weightVals.push_back(envelope); // weights could be sparse, fill in default weight if no data
				weightIndex++;
			}
			MDataHandle value = hWeights.inputValue(&status);
			weightVals.push_back(value.asFloat() * envelope);
			weightIndex++;
		}
		// now we have written the last non-zero weight into weightVals, but the last non-zero weight
		// doesn't have to be for the last vertex in the buffer. Add more default values if necessary.
		weightVals.resize(numVerts, envelope);
	}
	return MStatus::kSuccess;
}

MStatus BlurRelax::initialize() {
	MStatus status;
	MFnEnumAttribute eAttr;
	MFnNumericAttribute nAttr;

	aBorderBehavior = eAttr.create("borderBehavior", "bb", BB_SLIDE, &status);
	CHECKSTAT("aBorderBehavior");
	eAttr.setKeyable(false);
	eAttr.setChannelBox(true);
	status = eAttr.addField("None", BB_NONE);
	CHECKSTAT("aBorderBehavior");
	status = eAttr.addField("Pin", BB_PIN);
	CHECKSTAT("aBorderBehavior");
	status = eAttr.addField("Slide", BB_SLIDE);
	CHECKSTAT("aBorderBehavior");
	status = addAttribute(aBorderBehavior);
	CHECKSTAT("aBorderBehavior");
	status = attributeAffects(aBorderBehavior, outputGeom);
	CHECKSTAT("aBorderBehavior");

	aHardEdgeBehavior = eAttr.create("hardEdgeBehavior", "hb", HB_SLIDE, &status);
	CHECKSTAT("aHardEdgeBehavior");
	eAttr.setKeyable(false);
	eAttr.setChannelBox(true);
	status = eAttr.addField("None", HB_NONE);
	CHECKSTAT("aHardEdgeBehavior");
	status = eAttr.addField("Pin", HB_PIN);
	CHECKSTAT("aHardEdgeBehavior");
	status = eAttr.addField("Slide", HB_SLIDE);
	CHECKSTAT("aHardEdgeBehavior");
	status = addAttribute(aHardEdgeBehavior);
	CHECKSTAT("aHardEdgeBehavior");
	status = attributeAffects(aHardEdgeBehavior, outputGeom);
	CHECKSTAT("aHardEdgeBehavior");

	aGroupEdgeBehavior = eAttr.create("groupEdgeBehavior", "gb", GB_PIN, &status);
	CHECKSTAT("aGroupEdgeBehavior");
	eAttr.setKeyable(false);
	eAttr.setChannelBox(true);
	status = eAttr.addField("None", GB_NONE);
	CHECKSTAT("aGroupEdgeBehavior");
	status = eAttr.addField("Pin", GB_PIN);
	CHECKSTAT("aGroupEdgeBehavior");
	status = eAttr.addField("Slide", GB_SLIDE);
	CHECKSTAT("aGroupEdgeBehavior");
	status = addAttribute(aGroupEdgeBehavior);
	CHECKSTAT("aGroupEdgeBehavior");
	status = attributeAffects(aGroupEdgeBehavior, outputGeom);
	CHECKSTAT("aGroupEdgeBehavior");

	aReproject = nAttr.create("reproject", "rp", MFnNumericData::kBoolean, false, &status);
	CHECKSTAT("aReproject");
	nAttr.setKeyable(false);
	nAttr.setChannelBox(true);
	status = addAttribute(aReproject);
	CHECKSTAT("aReproject");
	status = attributeAffects(aReproject, outputGeom);
	CHECKSTAT("aReproject");

	// Taubin Bias is divided by 1000 internally
	aTaubinBias = nAttr.create("preserveVolume", "pv", MFnNumericData::kFloat, 0.0f, &status);
	CHECKSTAT("aTaubinBias");
	nAttr.setMin(0.0f);
	nAttr.setMax(2.0f);
	nAttr.setChannelBox(true);
	status = addAttribute(aTaubinBias);
	CHECKSTAT("aTaubinBias");
	status = attributeAffects(aTaubinBias, outputGeom);
	CHECKSTAT("aTaubinBias");

	aIterations = nAttr.create("iterations", "i", MFnNumericData::kFloat, 0, &status);
	CHECKSTAT("aIterations");
	nAttr.setMin(0.0);
	nAttr.setChannelBox(true);
	status = addAttribute(aIterations);
	CHECKSTAT("aIterations");
	status = attributeAffects(aIterations, outputGeom);
	CHECKSTAT("aIterations");

	MGlobal::executeCommand("makePaintable -attrType \"multiFloat\" -sm \"deformer\" \"" DEFORMER_NAME "\" \"weights\";");

	return MStatus::kSuccess;
}

MStatus BlurRelax::deform(MDataBlock& dataBlock, MItGeometry& vertIter, const MMatrix& matrix, UINT multiIndex) {
	// Short circuit if the envelope is zero
	MStatus status;
	MDataHandle hEnv = dataBlock.inputValue(envelope);
	float env = hEnv.asFloat();

	MDataHandle hIter = dataBlock.inputValue(aIterations);
	float iterations = hIter.asFloat();

	if (iterations > 0 && env > 0.0f){
		// Get the data from the node
		MDataHandle hBorder = dataBlock.inputValue(aBorderBehavior);
		short bb = hBorder.asShort();
		MDataHandle hHardEdge = dataBlock.inputValue(aHardEdgeBehavior);
		short hb = hHardEdge.asShort();
		MDataHandle hGroupEdge = dataBlock.inputValue(aGroupEdgeBehavior);
		short gb = hGroupEdge.asShort();
		MDataHandle hReproject = dataBlock.inputValue(aReproject);
		bool reproject = hReproject.asBool();

		MDataHandle hTBias = dataBlock.inputValue(aTaubinBias);
		float tBias = hTBias.asFloat();
		// volume preservation uses 2 steps per iteration
		// so half the number of iterations if I'm volume preserving
		// The iterations interpolating as floats takes care of 99% of the jumping
		// There *will* be a tiny jump from 0.0 on preserve volume if
		// reprojection is also turned on. I don't think I can get rid of that one
		if (tBias > 0.0f) {
			iterations /= 2.0f;
		}
		// So the numbers shown are intuitive to the user
		// 0 maps to 1.0 and 1 maps to -1.05
		// -1.05 is used because taubin smoothing needs something
		// just a little over 1.0 to truly preserve volume,
		// and that value looked good on my test mesh.
		tBias = -2.05f * tBias + 1.0f;

		// get the input mesh corresponding to this output
		MObject thisNode = this->thisMObject();
		MPlug inPlug(thisNode, input);
		inPlug.selectAncestorLogicalIndex(multiIndex, input);
		MDataHandle hInput = dataBlock.inputValue(inPlug);
		MObject mesh = hInput.asMesh();

		// Get the point values
		MFnMesh meshFn(mesh);
		UINT numVerts = meshFn.numVertices();

		// Get the group data
		std::vector<char> group;
		group.resize(numVerts);
		for (; !vertIter.isDone(); vertIter.next()) {
			group[vertIter.index()] = true;
		}
		vertIter.reset();

		// Get the true smooth-edge data
		std::vector<char> trueSmoothEdges;
		trueSmoothEdges.resize(meshFn.numEdges());
		MItMeshEdge edgeIter(mesh);
		for (; !edgeIter.isDone(); edgeIter.next()) {
			if (edgeIter.isSmooth()) {
				trueSmoothEdges[edgeIter.index()] = true;
			}
		}
		edgeIter.reset();

		// Hash the topology to see if we need to re-build the data
		// I could be smarter about this and only recompute things when needed
		// however, those re-computations wouldn't save frames when it mattered
		MIntArray mfaceCounts, mvertIdxs;
		int *faceCounts, *vertIdxs;
		meshFn.getVertices(mfaceCounts, mvertIdxs);
		faceCounts = new int[mfaceCounts.length()];
		vertIdxs = new int[mvertIdxs.length()];
		unsigned long long tFaceHash = XXH64(faceCounts, mfaceCounts.length()*sizeof(int), 0);
		unsigned long long tEdgeHash = XXH64(vertIdxs, mvertIdxs.length()*sizeof(int), 0);
		unsigned long long tGroupHash = XXH64(group.data(), group.size()*sizeof(char), 0);
		unsigned long long tSmoothHash = XXH64(trueSmoothEdges.data(), trueSmoothEdges.size()*sizeof(char), 0);
		delete[] faceCounts;
		delete[] vertIdxs;

		if (
			(bbCheck != bb) || (hbCheck != hb) || (gbCheck != gb) ||
			(tFaceHash != faceHash) || (tEdgeHash != edgeHash) ||
			(tGroupHash != groupHash) || (tSmoothHash != smoothHash)
		) {
			bbCheck = bb; hbCheck = hb; gbCheck = gb;
			faceHash = tFaceHash; edgeHash = tEdgeHash;
			groupHash = tGroupHash; smoothHash = tSmoothHash;

			// Populate the variables with *SPECIALLY ORDERED* data
			// all vertex data is now shuffled by the order vector
			buildQuickData(mesh, vertIter, bb, hb, gb, reproject, group);
		}

		// This can happen if the user is pinning all the points
		// or all the edges are hard (like when you import an obj)
		if (neighbors.empty()) {
			return status;
		}
		
		// Build the raw float data buffers
		pointArray_t mpa;
		float_t(*reoVerts)[4] = new float_t[numVerts][4];
		meshFn.getPoints(mpa);

		for (size_t i = 0; i < numVerts; ++i) {
			reoVerts[i][0] = mpa[order[i]].x;
			reoVerts[i][1] = mpa[order[i]].y;
			reoVerts[i][2] = mpa[order[i]].z;
			reoVerts[i][3] = mpa[order[i]].w;
		}

		// Calculate the relax, and store in verts
		quickRelax(mesh, bb, hb, gb, reproject, tBias, iterations, numVerts, group, reoVerts);

		// Get the painted weight values
		std::vector<float> weightVals;
		status = getTrueWeights(mesh, dataBlock, multiIndex, weightVals, env);

		// Finally set the output
		for (; !vertIter.isDone(); vertIter.next()) {
			const UINT idx = vertIter.index();
			vertIter.setPosition((weightVals[idx]) * point_t(reoVerts[invOrder[idx]]) + (1.0 - weightVals[idx]) * vertIter.position());
		}

		// Make sure to clean up after myself
		delete [] reoVerts;
	}
	return status;
}

MStatus initializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj, "Blur Studio", "1.0", "Any");
	result = plugin.registerNode(DEFORMER_NAME, BlurRelax::id, BlurRelax::creator, BlurRelax::initialize, MPxNode::kDeformerNode);

	return result;
}

MStatus uninitializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj);
	result = plugin.deregisterNode(BlurRelax::id);
	return result;
}

void BlurRelax::buildQuickData(
	MObject &mesh,
	MItGeometry& vertIter,
	short borderBehavior,
	short hardEdgeBehavior,
	short groupEdgeBehavior,
	bool reproject,
	std::vector<char> &group
) {
	// This takes the mesh, Gets all the required data, reorders so
	// the verts are in descending order of valence, and gets the subgroup
	// I'm trying to pre-process everything I can right here so I don't have to
	// branch in quickLaplacianSmooth() so auto-vectorization works
	MFnMesh meshFn(mesh);
	UINT numVertices = meshFn.numVertices();

	std::vector<char> rawPinPoints;
	std::vector<UINT> rawCreaseCount;
	rawCreaseCount.resize(numVertices);
	rawPinPoints.resize(numVertices);

	std::vector<UINT> groupBorders;
	std::vector<char> isGroupBorder;
	isGroupBorder.resize(numVertices);

	std::vector<std::vector<UINT>> rawNeighbors;
	std::vector<std::vector<char>> rawHardEdges;
	std::vector<std::vector<char>> rawPinBorders;
	rawHardEdges.resize(numVertices);
	rawNeighbors.resize(numVertices);

	// Get the connectivity data

	// Instead of treating neighbors as an un-directed graph, I treat it as a directed graph
	// This way I can remove influences from vertices, but not their neighbors
	// So I can remove all neighbors from a vertex to pin it, but its neighbors can still use it as an influence.
	// The big trick is that all that pinning computation only happens once, right here.
	// And after all this junk, I only have to apply one set of very simple, very optimizable rules

	MItMeshEdge edgeIter(mesh);
	for (; !edgeIter.isDone(); edgeIter.next()) {
		const UINT start = edgeIter.index(0);
		const UINT end = edgeIter.index(1);
		// Hard means the edge has a special behavior applied to it
		// Pinned means that special behavior is to not move
		// hard and unpinned means slide if there are exactly 2 connected hard edges
		// hard and unpinned means pin otherwise
		// TODO: Also handle *IGNORED* neighbors
		const bool onBound = edgeIter.onBoundary();
		const bool isHard = !edgeIter.isSmooth();
		const bool hard = (onBound && borderBehavior) || (isHard && !onBound && hardEdgeBehavior);
		const bool pin = ((borderBehavior == BB_PIN) && onBound) || ((hardEdgeBehavior == HB_PIN) && isHard && !onBound);

		const bool startInGroup = group[start];
		const bool endInGroup = group[end];
		if (!startInGroup && !endInGroup) continue;

		if (startInGroup && endInGroup) {
			if (pin) {
				rawPinPoints[start] = true;
				rawPinPoints[end] = true;
			}
			else {
				rawNeighbors[start].push_back(end);
				rawHardEdges[start].push_back(hard);
				rawNeighbors[end].push_back(start);
				rawHardEdges[end].push_back(hard);
				if (hard) {
					++rawCreaseCount[start];
					++rawCreaseCount[end];
				}
			}
		}
		else if (!startInGroup) {
			groupBorders.push_back(end);
			isGroupBorder[end] = true;
		}
		else if (!endInGroup) {
			groupBorders.push_back(start);
			isGroupBorder[start] = true;
		}
	}

	const bool gbHard = groupEdgeBehavior;
	const bool gbPin = groupEdgeBehavior == GB_PIN;
	for (UINT gb: groupBorders){
		const auto nnn = rawNeighbors[gb];
		for (size_t i = 0; i < nnn.size(); ++i) {
			const UINT n = nnn[i];
			if (isGroupBorder[n]  ){
				rawHardEdges[gb][i] = rawHardEdges[gb][i] || gbHard;
				if (gbPin){
					rawPinPoints[gb] = true;
					rawPinPoints[n] = true;
				}
				else if (gbHard && (gb > n)){
					// Each edge this way is visited twice
					// So check (gb>n) to only increment once
					++rawCreaseCount[gb];
					++rawCreaseCount[n];
				}
			}
		}
	}

	// if a vert has hard neighbors, remove all non-hard neighbors
	// if a vert is pinned remove all of its neighbors
	// if a vert has a pinned neighbor, remove it

	std::vector<float_t> rawShiftVal;
	rawShiftVal.resize(numVertices);
	std::fill(rawShiftVal.begin(), rawShiftVal.end(), 0.5);

	for (size_t i = 0; i < rawNeighbors.size(); ++i) {
		if ((rawCreaseCount[i] != 0) || rawPinPoints[i]) {
			if (rawCreaseCount[i] != 2)
				rawPinPoints[i] = true;

			std::vector<UINT> newNeigh;
			std::vector<char> newHard;
			if (!rawPinPoints[i]) {
				for (size_t j = 0; j < rawNeighbors[i].size(); ++j) {
					if (rawHardEdges[i][j]) {
						newNeigh.push_back(rawNeighbors[i][j]);
						newHard.push_back(rawHardEdges[i][j]);
					}
				}
			}
			rawNeighbors[i] = newNeigh;
			rawHardEdges[i] = newHard;
			rawShiftVal[i] = 0.25;
		}
	}

	// Sort the verts by their valence count
	order.resize(numVertices);
	std::iota(order.begin(), order.end(), 0u); // like python range()

	// compare using > for a sort from high to low
	std::sort(order.begin(), order.end(), [&rawNeighbors](size_t a, size_t b) {return rawNeighbors[a].size() > rawNeighbors[b].size(); });

	// Build the "transposed" neighbor and hard edge values
	size_t maxValence = rawNeighbors[order[0]].size();
	neighbors.clear();
	neighbors.resize(maxValence);
	hardEdges.clear();
	hardEdges.resize(maxValence);
	creaseCount.resize(numVertices);
	pinPoints.resize(numVertices);

	valence.resize(numVertices*4);
	shiftVal.resize(numVertices*4);

	invOrder.resize(order.size());

	for (size_t i = 0; i < order.size(); ++i) {
		invOrder[order[i]] = i;
	}

	for (size_t i = 0; i < numVertices; ++i) {
		const auto &neigh = rawNeighbors[order[i]];
		const auto &hards = rawHardEdges[order[i]];
		size_t vale = neigh.size();
		for (size_t n = 0; n < vale; ++n) {
			neighbors[n].push_back(invOrder[neigh[n]]);
			hardEdges[n].push_back(invOrder[hards[n]]);
		}
		creaseCount[i] = rawCreaseCount[order[i]];
		pinPoints[i] = rawPinPoints[order[i]];

		// Vectorizing flattens the vert list, so
		// I need these per vert, per component
		for (size_t xx = 0; xx < 4; ++xx) {
			valence[4 * i + xx] = float_t(vale);
			shiftVal[4 * i + xx] = rawShiftVal[order[i]];
		}
	}
}

void BlurRelax::quickRelax(
	MObject &mesh,
	const short borderBehavior,
	const short hardEdgeBehavior,
	const short groupEdgeBehavior,
	const bool reproject,
	const float taubinBias,
	const float_t iterations,
	const UINT numVerts,
	const std::vector<char> &group,
	float_t(*verts)[4]
) {
	bool rpEdges = (borderBehavior == BB_SLIDE) || (hardEdgeBehavior == HB_SLIDE) || (groupEdgeBehavior == GB_SLIDE);
	std::vector<size_t> groupIdxs;

	float_t (*baseVerts)[4];
	if (rpEdges) {
		for (size_t i = 0; i < group.size(); ++i) {
			if (group[i]) groupIdxs.push_back(i);
		}
		// make a copy of the original verts only if they'll be used for edge reprojection
		baseVerts = new float_t[numVerts][4];
		memcpy(&(baseVerts[0][0]), &(verts[0][0]), 4 * numVerts * sizeof(float_t));
	}

	float_t(*prevVerts)[4];
	prevVerts = new float_t[numVerts][4];

	float_t iterT, iterFI;
	iterT = modf(iterations, &iterFI);
	UINT iterI = (UINT)iterFI;
	if (iterT > 0.0) {
		iterI += 1;
	}
	
	size_t nonzeroValence = neighbors[0].size();

	MStatus status;
	MMeshIntersector octree;
	MObject smoothMeshPar, smoothMesh;
	MFnMesh meshFn(mesh);
	if (reproject) {
		MFnMeshData smoothMeshParFn;
		MMeshSmoothOptions smoothOpt;
		smoothMeshPar = smoothMeshParFn.create();
		smoothOpt.setDivisions(1);
		smoothOpt.setKeepBorderEdge(rpEdges);
		smoothOpt.setSubdivisionType(MMeshSmoothOptions::kCatmullClark);
		smoothMesh = meshFn.generateSmoothMesh(smoothMeshPar, &smoothOpt);
		octree.create(smoothMesh);
	}

	for (size_t r = 0; r < iterI; ++r) {
		if ((r == iterI - 1) && (iterT > 0.0)){
			// Store the next-to-last iteration to interpolate with
			memcpy(&(prevVerts[0][0]), &(verts[0][0]), 4 * numVerts * sizeof(float_t));
		}
		quickLaplacianSmooth(verts, numVerts, neighbors, valence, shiftVal, pinPoints);
		if (taubinBias < 1.0){
			quickLaplacianSmooth(verts, numVerts, neighbors, valence, shiftVal, pinPoints, taubinBias);
		}

		if (rpEdges) {
			edgeProject(baseVerts, groupIdxs, invOrder, neighbors, hardEdges, creaseCount, verts);
		}

		if (reproject) {
			#pragma omp parallel for if(numVerts>2000)
			for (int i = 0; i < nonzeroValence; ++i) {
				if ((creaseCount[i] == 0) && (group[order[i]])) {
					point_t mf(verts[i][0], verts[i][1], verts[i][2]);
					MPointOnMesh pom;
					octree.getClosestPoint(mf, pom);
					point_t gpf = pom.getPoint();
					//mfp.set(gpf, i);
					verts[i][0] = gpf[0];
					verts[i][1] = gpf[1];
					verts[i][2] = gpf[2];
				}
			}
		}
	}

	// Interpolate between prevVerts and verts based on iterT
	if (iterT > 0.0) {
		// This should vectorize
		float_t * vv = &verts[0][0];
		float_t * pv = &prevVerts[0][0];
		for (size_t i = 0; i < numVerts * 4; ++i) {
			vv[i] = ((vv[i] - pv[i]) * iterT) + pv[i];
		}
	}

	if (rpEdges) delete[] baseVerts;
	delete[] prevVerts;
}

