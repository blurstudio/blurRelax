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

#define CHECKSTAT(m) if (!status) {status.perror(m); return status;};

#define BB_NONE   0
#define BB_PIN    1
#define BB_SLIDE  2

#define HB_NONE   0
#define HB_PIN    1
#define HB_SLIDE  2

void edgeProject(
	const float basePoints[][4],
	const std::vector<size_t> &group,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<std::vector<bool>> &hardEdges,
	const std::vector<UINT> &creaseCount,
	float smoothPoints[][4]
) {
	/*
	This code takes the edges that have been defined as hard
	and reprojects them onto the original neighboring hard edges

	Perhaps in the future, I can project onto non-branching
	strings of edges. But for now, it's just direct neighbors
	*/

	std::vector<UINT> neigh;

	for (size_t gidx = 0; gidx < group.size(); ++gidx) {
		// If we have "hard edges" we have already removed
		// any soft edges, so if the crease count for this vertex is zero
		// we can just completely skip it. And if it's >0, then we can
		// assume that all stored neighbors are hard.
		size_t idx = group[gidx];

		float *avg = smoothPoints[idx];
		const float *basePos = basePoints[idx];
		if (creaseCount[idx] != 2) {
			//for (size_t x = 0; x < 3; ++x)
				//avg[x] = basePos[x];
			continue;
		}
		float keep[3], delta[3], edge[3];

		neigh.clear();
		for (size_t i = 0; i < neighbors.size(); ++i) {
			if (idx < neighbors[i].size()) {
				neigh.push_back(neighbors[i][idx]);
			}
			else break;
		}

		// TODO: Find the hard edge strings and reproject onto those
		float minLen = std::numeric_limits<float>::max();
		bool found = false;

		// avg - basePos
		for (size_t x = 0; x < 3; ++x)
			delta[x] = avg[x] - basePos[x];

		for (size_t i = 0; i < neigh.size(); ++i) {
			UINT n = neigh[i];

			// normalized(prevPoints[n] - basePos)
			for (size_t x = 0; x < 3; ++x)
				edge[x] = basePoints[n][x] - basePos[x];
			float elen = sqrt(edge[0] * edge[0] + edge[1] * edge[1] + edge[2] * edge[2]);
			for (size_t x = 0; x < 3; ++x)
				edge[x] /= elen;

			// dot(delta, edge)
			float dd = delta[0] * edge[0] + delta[1] * edge[1] + delta[2] * edge[2];

			float dn[3] = { 0.0f, 0.0f, 0.0f };
			if (dd > 0.0) {
				for (size_t x = 0; x < 3; ++x)
					dn[x] = edge[x] * dd;
			}
			float xx[3];

			// delta - dn
			for (size_t x = 0; x < 3; ++x)
				xx[x] = delta[x] - dn[x];

			// dot(xx, xx)
			float len2 = xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2];

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
	float verts[][4],
	const size_t numVerts,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<float> &valence,
	const std::vector<float> &shiftVal,
	const std::vector<float> &shiftComp,
	const std::vector<bool> &pinPoints
) {
	/*
	All the crazy hoops I've jumped through are to make auto-vectorization work
	
	neighbors and neighborOffsets work together to be a kind of pseudo-transpose
	of a std::vector of neighbors per vert with that vector sorted so the verts
	with the most neighbors were at the top.
	neighborOffsets contains the index offsets of vertices with at least [index] of neighbors
	So a single-subdivided cube would have 8 3-valence, and 18 4-valence for a total of 26 verts, and 96 neighbors
	So the neighborOffsets would be [0, 26, 52, 78, 96] meaning that
	verts 0-25 have at least 1 neighbor
	verts 26-51 have at least 2 neighbors
	verts 52-77 have at least 3 neighbors
	verts 78-95 have at least 4 neighbors
	*/

	// The __restrict keyword tells the compiler that *output and all it's offsets
	// are not pointed to by any other pointer in this scope
	// This allows for auto-vectorization
	// We assume that outpuoutput has been resized properly

	float * __restrict outComp = new float[numVerts];
	float * __restrict comp = new float[numVerts];
	float(*output)[4] = new float[numVerts][4];

	// number of nonzero valence
	size_t nzv = neighbors[0].size();

	// loop once for each component x,y,z
	for (size_t c = 0; c < 3; ++c) {
		// copy the component into contiguous memory
		for (size_t i = 0; i < numVerts; ++i) {
			comp[i] = verts[i][c];
		}
		// zero the accumulator
		for (size_t i = 0; i < numVerts; ++i) {
			outComp[i] = 0.0f;
		}

		for (size_t ncIdx = 0; ncIdx < neighbors.size(); ++ncIdx) {
			const auto &nCol = neighbors[ncIdx];
			size_t nColCount = nCol.size();
			for (size_t i = 0; i < nColCount; ++i) {
				outComp[i] = outComp[i] + comp[nCol[i]];
			}
		}

		// Depending on the compiler optimization, it may be faster to break up this line
		// Gotta test
		for (size_t i = 0; i < nzv; ++i) {
			outComp[i] = shiftVal[i] * (outComp[i] / valence[i]) + shiftComp[i] * comp[i];
		}

		// reset the data into the output
		for (size_t i = 0; i < numVerts; ++i) {
			output[i][c] = outComp[i];
		}
	}

	for (size_t i = 0; i < nzv; ++i) {
		//if (!pinPoints[i]) {
			for (size_t x = 0; x < 3; ++x)
				verts[i][x] = output[i][x];
		//}
	}

	delete outComp;
	delete comp;
	delete output;
}

class BlurRelax : public MPxDeformerNode {
	public:
		BlurRelax();
		~BlurRelax() override;

		static void* creator();
		static MStatus initialize();
		MStatus compute(const MPlug& plug, MDataBlock& dataBlock) override;

	public:
		// local node attributes
		static MObject aIterations;
		static MObject aBorderBehavior;
		static MObject aHardEdgeBehavior;
		static MObject aPinGroupEdges;
		static MObject aReproject;

		static MTypeId id;
	private:

		MStatus getTrueWeights(
				MObject &mesh,
				MDataBlock &dataBlock,
				UINT index,
				std::vector<float> &weightVals,
				float envelope
				) const;

		void BlurRelax::buildQuickData(
			MDataHandle &meshHandle,
			UINT groupId,
			short borderBehavior,
			short hardEdgeBehavior,
			bool reproject,
			bool pinGroupEdges,
			std::vector<bool> &group,
			std::vector<size_t> &order,
			std::vector<std::vector<size_t>> &neighbors,
			std::vector<std::vector<bool>> &hardEdges,
			std::vector<float> &shiftVal, // normally 0.5, but it's 0.25 if on a hard edge
			std::vector<float> &shiftComp, // normally 0.5, but it's 0.75 if on a hard edge
			std::vector<float> &valence, // as float for vectorizing
			std::vector<bool> &pinPoints,
			std::vector<UINT> &creaseCount,
			float(*verts)[4]
		);

		void quickRelax(
			const MObject &mesh,
			const short borderBehavior,
			const short hardEdgeBehavior,
			const bool reproject,
			const bool pinGroupEdges,
			const UINT iterations,
			const UINT numVerts,
			const std::vector<bool> &group,
			const std::vector<size_t> &order,
			const std::vector<std::vector<size_t>> &neighbors,
			const std::vector<std::vector<bool>> &hardEdges,
			const std::vector<float> &shiftVal, // normally 0.5, but it's 0.25 if on a hard edge
			const std::vector<float> &shiftComp, // normally 0.5, but it's 0.75 if on a hard edge
			const std::vector<float> &valence, // as float for vectorizing
			const std::vector<bool> &pinPoints,
			const std::vector<UINT> &creaseCount,
			float(*verts)[4], // already resized
			const float(*baseVerts)[4] // already resized
		);
};

MTypeId BlurRelax::id(0x001226FD);

// local attributes
MObject BlurRelax::aIterations;
MObject BlurRelax::aBorderBehavior;
MObject BlurRelax::aHardEdgeBehavior;
MObject BlurRelax::aPinGroupEdges;
MObject BlurRelax::aReproject;

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
	// it is possible that the jumpToElement fails.  In that case all weights are 1.
	if (!status) {
		for (unsigned int i = 0; i < meshFn.numVertices(); i++)
			weightVals.push_back(envelope);
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
				weightVals.push_back(0.0f); // weights could be sparse, fill in zero weight if no data
				weightIndex++;
			}
			MDataHandle value = hWeights.inputValue(&status);
			weightVals.push_back(value.asFloat() * envelope);
			weightIndex++;
		}
		// now we have written the last non-zero weight into weightVals, but the last non-zero weight
		// doesn't have to be for the last vertex in the buffer.  Add more zero values if necessary.
		while (weightIndex < numVerts) {
			weightVals.push_back(0.0f); // weights could be sparse, fill in zero weight if no data
			weightIndex++;
		}
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

	aPinGroupEdges = nAttr.create("pinGroupEdges", "pge", MFnNumericData::kBoolean, true, &status);
	CHECKSTAT("aPinGroupEdges");
	nAttr.setKeyable(false);
	nAttr.setChannelBox(true);
	status = addAttribute(aPinGroupEdges);
	CHECKSTAT("aPinGroupEdges");
	status = attributeAffects(aPinGroupEdges, outputGeom);
	CHECKSTAT("aPinGroupEdges");

	aReproject = nAttr.create("reproject", "rp", MFnNumericData::kBoolean, true, &status);
	CHECKSTAT("aReproject");
	nAttr.setKeyable(false);
	nAttr.setChannelBox(true);
	status = addAttribute(aReproject);
	CHECKSTAT("aReproject");
	status = attributeAffects(aReproject, outputGeom);
	CHECKSTAT("aReproject");

	aIterations = nAttr.create("iterations", "i", MFnNumericData::kInt, 10, &status);
	CHECKSTAT("aIterations");
	nAttr.setMin(0);
	nAttr.setChannelBox(true);
	status = addAttribute(aIterations);
	CHECKSTAT("aIterations");
	status = attributeAffects(aIterations, outputGeom);
	CHECKSTAT("aIterations");

	return MStatus::kSuccess;
}

MStatus BlurRelax::compute(const MPlug& plug, MDataBlock& dataBlock) {
	MStatus status = MStatus::kUnknownParameter;
	if (plug.attribute() == outputGeom) {
		// Get the behavior options
		MDataHandle hBorder = dataBlock.inputValue(aBorderBehavior);
		short bb = hBorder.asShort();
		MDataHandle hHardEdge = dataBlock.inputValue(aHardEdgeBehavior);
		short hb = hHardEdge.asShort();
		MDataHandle hPGE = dataBlock.inputValue(aPinGroupEdges);
		bool pge = hPGE.asBool();
		MDataHandle hReproject = dataBlock.inputValue(aReproject);
		bool reproject = hReproject.asBool();
		MDataHandle hIter = dataBlock.inputValue(aIterations);
		int iterations = hIter.asInt();
		MDataHandle hEnv = dataBlock.inputValue(envelope);
		float env = hEnv.asFloat();

		// get the input corresponding to this output
		unsigned int index = plug.logicalIndex();
		MObject thisNode = this->thisMObject();
		MPlug inPlug(thisNode, input);
		inPlug.selectAncestorLogicalIndex(index, input);
		MDataHandle hInput = dataBlock.inputValue(inPlug);

		// get the input geometry and input groupId
		MDataHandle hGeom = hInput.child(inputGeom);
		MDataHandle hGroup = hInput.child(groupId);
		unsigned int groupId = hGroup.asLong();

		MDataHandle hOutput = dataBlock.outputValue(plug);
		hOutput.copy(hGeom);

		// Maya might do some tricky stuff like not store the weight array at all for certain weight
		// values so we can't count on an array existing in the weightList.  For the OpenCL Kernel
		// we want an array with one weight in it per vertex, we need to build it carefully here.

		// Two possibilities: we could have a sparse array in weightList[multiIndex] or there could be nothing in weightList[multiIndex].
		// if nothing is there then all the weights at 1.0f.

		// Get a handle to the weight array we want.
		if (env > 0.0f){
			MObject mesh = hOutput.asMesh();
			MFnMesh meshFn(mesh);
			MFloatPointArray mfp;
			meshFn.getPoints(mfp);
			UINT numVerts = meshFn.numVertices();

			float (*verts)[4] = new float[numVerts][4];
			float (*baseVerts)[4] = new float[numVerts][4];
			mfp.get(verts);
			
			std::vector<float> weightVals;
			status = getTrueWeights(mesh, dataBlock, index, weightVals, env);

			//shiftVal, // normally 0.5, but it's 0.25 if on a hard edge
			//shiftComp, // normally 0.5, but it's 0.75 if on a hard edge
			//valence, // as float for vectorizing

			std::vector<bool> group;
			std::vector<size_t> order;
			std::vector<std::vector<size_t>> neighbors;
			std::vector<std::vector<bool>> hardEdges;
			std::vector<float> shiftVal; // normally 0.5; but it's 0.25 if on a hard edge
			std::vector<float> shiftComp; // normally 0.5; but it's 0.75 if on a hard edge
			std::vector<float> valence; // as float for vectorizing
			std::vector<bool> pinPoints;
			std::vector<UINT> creaseCount;

			buildQuickData(hOutput, groupId, bb, hb, reproject, pge,
				group, order, neighbors, hardEdges, shiftVal, shiftComp, valence,
				pinPoints, creaseCount, verts);

			for (size_t i = 0; i < numVerts; ++i) {
				for (size_t x = 0; x < 3; ++x) baseVerts[i][x] = verts[i][x];
			}

			quickRelax(mesh, bb, hb, reproject, pge,
				iterations, numVerts, group, order, neighbors, hardEdges, shiftVal,
				shiftComp, valence, pinPoints, creaseCount, verts, baseVerts);

			// undo ordering
			for (size_t i = 0; i < order.size(); ++i) {
				for (size_t x = 0; x < 3; ++x)
					baseVerts[order[i]][x] = verts[i][x];
				baseVerts[order[i]][3] = 1.0f; // maya keeps 
			}

			MFloatPointArray outMfp(baseVerts, numVerts);
			meshFn.setPoints(outMfp);

			// TODO: vectorize this??
			MItGeometry vertIter(hOutput, groupId, false);
			for (; !vertIter.isDone(); vertIter.next()) {
				UINT idx = vertIter.index();
				vertIter.setPosition((1.0 - weightVals[idx]) * mfp[idx] + weightVals[idx] * vertIter.position());
			}

			delete verts;
			delete baseVerts;
		}
		hOutput.setClean();
	}
	return status;
}

MStatus initializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj, "Blur Studio", "1.0", "Any");
	result = plugin.registerNode("BlurRelax", BlurRelax::id, BlurRelax::creator, BlurRelax::initialize, MPxNode::kDeformerNode);

	return result;
}

MStatus uninitializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj);
	result = plugin.deregisterNode(BlurRelax::id);
	return result;
}

void BlurRelax::buildQuickData(
	MDataHandle &meshHandle,
	UINT groupId,
	short borderBehavior,
	short hardEdgeBehavior,
	bool reproject,
	bool pinGroupEdges,
	std::vector<bool> &group,
	std::vector<size_t> &order,
	std::vector<std::vector<size_t>> &neighbors,
	std::vector<std::vector<bool>> &hardEdges,
	std::vector<float> &shiftVal, // normally 0.5, but it's 0.25 if on a hard edge
	std::vector<float> &shiftComp, // normally 0.5, but it's 0.75 if on a hard edge
	std::vector<float> &valence, // as float for vectorizing
	std::vector<bool> &pinPoints,
	std::vector<UINT> &creaseCount,
	float(*verts)[4] // already resized
) {
	// This takes the mesh, Gets all the required data, reorders so
	// the verts are in descending order of valence, and gets the subgroup
	// I'm trying to pre-process everything I can right here so I don't have to
	// branch in quickLaplacianSmooth() so auto-vectorization works

	MObject mesh = meshHandle.asMesh();
	MFnMesh meshFn(mesh);
	UINT numVertices = meshFn.numVertices();
	float(*rawVerts)[4] = new float[numVertices][4];

	for (size_t i = 0; i < numVertices; ++i) {
		for (size_t x = 0; x < 3; ++x)
			rawVerts[i][x] = verts[i][x];
	}

	std::vector<bool> rawPinPoints;
	std::vector<UINT> rawCreaseCount;
	rawCreaseCount.resize(numVertices);
	rawPinPoints.resize(numVertices);

	std::vector<std::vector<UINT>> rawNeighbors;
	std::vector<std::vector<bool>> rawHardEdges;
	std::vector<std::vector<bool>> rawPinBorders;
	rawHardEdges.resize(numVertices);
	rawNeighbors.resize(numVertices);

	// Get the group data
	group.resize(numVertices);
	MItGeometry vertIter(meshHandle, groupId, true);
	for (; !vertIter.isDone(); vertIter.next()) {
		group[vertIter.index()] = true;
	}

	// Get the connectivity data
	// The big trick is that all pinning happens once, right here
	// If a vertex has no neighbors, then it is never smoothed
	// so any vertices that are pinned just aren't given any neighbors
	// (Meaning if A has neighbor B, then B doesn't necessarily have neighbor A)
	MItMeshEdge edgeIter(mesh);
	for (; !edgeIter.isDone(); edgeIter.next()) {
		const UINT start = edgeIter.index(0);
		const UINT end = edgeIter.index(1);
		// Hard means the edge has a special behavior applied to it
		// Pinned means that special behavior is to not move
		// hard and unpinned means slide if there are exactly 2 connected hard edges
		// hard and unpinned means pin otherwise
		// TODO: Also handle *IGNORED* neighbors
		const bool hard = (edgeIter.onBoundary() && borderBehavior) || (!edgeIter.isSmooth() && !edgeIter.onBoundary() && hardEdgeBehavior);
		const bool pin = (borderBehavior == BB_PIN && edgeIter.onBoundary()) ||
			(hardEdgeBehavior == HB_PIN && !edgeIter.isSmooth() && !edgeIter.onBoundary());

		const bool startInGroup = group[start];
		const bool endInGroup = group[end];
		if (!startInGroup && !endInGroup) continue;

		if (startInGroup && endInGroup) {
			rawNeighbors[start].push_back(end);
			rawHardEdges[start].push_back(hard);
			rawNeighbors[end].push_back(start);
			rawHardEdges[end].push_back(hard);
			if (hard) {
				++rawCreaseCount[start];
				++rawCreaseCount[end];
			}
			if (pin) {
				rawPinPoints[start] = true;
				rawPinPoints[end] = true;
			}
		}
		else if (!startInGroup) {
			rawNeighbors[end].push_back(start);
			rawHardEdges[end].push_back(hard);
			if (hard) {
				++rawCreaseCount[start];
				++rawCreaseCount[end];
			}
			if (pinGroupEdges)
				rawPinPoints[start] = true;
		}
		else if (!endInGroup) {
			rawNeighbors[start].push_back(end);
			rawHardEdges[start].push_back(hard);
			if (hard) {
				++rawCreaseCount[end];
				++rawCreaseCount[start];
			}
			if (pinGroupEdges)
				rawPinPoints[end] = true;
		}
	}

	// if a vert has hard neighbors, remove all non-hard neighbors
	// if a vert is pinned remove all of its neighbors
	// if a vert has a pinned neighbor, remove it

	std::vector<float> rawShiftVal;
	std::vector<float> rawShiftComp;
	rawShiftVal.resize(numVertices);
	rawShiftComp.resize(numVertices);
	std::fill(rawShiftVal.begin(), rawShiftVal.end(), 0.5f);
	std::fill(rawShiftComp.begin(), rawShiftComp.end(), 0.5f);

	for (size_t i = 0; i < rawNeighbors.size(); ++i) {
		if (rawCreaseCount[i] != 0) {
			if (rawCreaseCount[i] != 2) rawPinPoints[i] = true;

			std::vector<UINT> newNeigh;
			std::vector<bool> newHard;
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
			rawShiftComp[i] = 0.75;
		}
	}

	// Sort the verts by their valence count
	order.resize(numVertices);
	std::iota(order.begin(), order.end(), 0u); // like python range()
	// compare using > for a sort from high to low
	std::sort(order.begin(), order.end(), [&rawNeighbors](UINT a, UINT b) {return rawNeighbors[a].size() > rawNeighbors[b].size(); });

	// Build the "transposed" neighbor and hard edge values
	valence.resize(numVertices);
	size_t maxValence = rawNeighbors[order[0]].size();
	neighbors.resize(maxValence);
	hardEdges.resize(maxValence);
	creaseCount.resize(numVertices);
	pinPoints.resize(numVertices);
	shiftVal.resize(numVertices);
	shiftComp.resize(numVertices);

	std::vector<size_t> invOrder(order.size());
	for (size_t i = 0; i < order.size(); ++i) {
		invOrder[order[i]] = i;
	}

	for (size_t i = 0; i < numVertices; ++i) {
		const std::vector<UINT> &neigh = rawNeighbors[order[i]];
		const std::vector<bool> &hards = rawHardEdges[order[i]];
		size_t vale = neigh.size();
		valence[i] = (float)vale;
		for (size_t n = 0; n < vale; ++n) {
			neighbors[n].push_back(invOrder[neigh[n]]);
			hardEdges[n].push_back(invOrder[hards[n]]);
		}
		for (size_t x = 0; x < 3; ++x)
			verts[i][x] = rawVerts[order[i]][x];
		creaseCount[i] = rawCreaseCount[order[i]];
		pinPoints[i] = rawPinPoints[order[i]];
		shiftVal[i] = rawShiftVal[order[i]];
		shiftComp[i] = rawShiftComp[order[i]];
	}
	delete rawVerts;
}

void BlurRelax::quickRelax(
	const MObject &mesh,
	const short borderBehavior,
	const short hardEdgeBehavior,
	const bool reproject,
	const bool pinGroupEdges,
	const UINT iterations,
	const UINT numVerts,
	const std::vector<bool> &group,
	const std::vector<size_t> &order,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<std::vector<bool>> &hardEdges,
	const std::vector<float> &shiftVal, // normally 0.5, but it's 0.25 if on a hard edge
	const std::vector<float> &shiftComp, // normally 0.5, but it's 0.75 if on a hard edge
	const std::vector<float> &valence, // as float for vectorizing
	const std::vector<bool> &pinPoints,
	const std::vector<UINT> &creaseCount,
	float(*verts)[4], // already resized
	const float(*baseVerts)[4] // already resized
) {

	std::vector<size_t> groupIdxs;
	for (size_t i = 0; i < group.size(); ++i) {
		if (group[i]) groupIdxs.push_back(i);
	}
	bool rpEdges = (borderBehavior == BB_SLIDE) || (hardEdgeBehavior == HB_SLIDE);

	size_t nonzeroValence = neighbors[0].size();


	MStatus status;
	MFnMesh meshFn(mesh);
	MMeshIntersector octree;
	MObject smoothMeshPar, smoothMesh;
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

	for (size_t r = 0; r < iterations; ++r) {
		quickLaplacianSmooth(verts, numVerts, neighbors, valence, shiftVal, shiftComp, pinPoints);
		if (rpEdges) {
			edgeProject(baseVerts, groupIdxs, neighbors, hardEdges, creaseCount, verts);
		}
		if (reproject) {
			#pragma omp parallel for if(numVerts>2000)
			for (int i = 0; i < nonzeroValence; ++i) {
				if ((creaseCount[i] == 0) && (group[order[i]])) {
					MFloatPoint mf(verts[i][0], verts[i][1], verts[i][2]);
					MPointOnMesh pom;
					octree.getClosestPoint(mf, pom);
					MFloatPoint gpf = pom.getPoint();
					//mfp.set(gpf, i);
					verts[i][0] = gpf[0];
					verts[i][1] = gpf[1];
					verts[i][2] = gpf[2];
				}
			}
		}
	}
}

