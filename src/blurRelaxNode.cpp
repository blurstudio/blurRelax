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
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MArrayDataHandle.h>
#include <maya/MGlobal.h>

#include <maya/MMeshIntersector.h>

#include <maya/MItMeshEdge.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItGeometry.h>
#include <maya/MPxDeformerNode.h> 

#include <maya/MFnEnumAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnComponent.h>
#include <maya/MFnGeometryFilter.h>

#include <maya/MPoint.h>
#include <maya/MFloatPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatVector.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MVector.h>
#include <maya/MVectorArray.h>

#include <maya/MMatrix.h>
#include <maya/MMatrixArray.h>

#include "blurRelaxNode.h"
#include "fastRelax.h"


MTypeId BlurRelax::id(0x001226FD);

// local attributes
MObject BlurRelax::aIterations;
MObject BlurRelax::aBorderBehavior;
MObject BlurRelax::aHardEdgeBehavior;
MObject BlurRelax::aGroupEdgeBehavior;
MObject BlurRelax::aReproject;
MObject BlurRelax::aReprojectDivs;
MObject BlurRelax::aTaubinBias;
MObject BlurRelax::aRecomputeTopo;
MObject BlurRelax::aDeltaMush;
MObject BlurRelax::aDeltaMult;
MObject BlurRelax::aDeltaBase;


BlurRelax::BlurRelax() {}
BlurRelax::~BlurRelax() {}
void* BlurRelax::creator() { return new BlurRelax(); }


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
	MFnTypedAttribute tAttr;

	aBorderBehavior = eAttr.create("borderBehavior", "bb", (UCHAR)B::SLIDE, &status);
	CHECKSTAT("aBorderBehavior");
	eAttr.setKeyable(false);
	eAttr.setChannelBox(true);
	status = eAttr.addField("None", (UCHAR)B::NONE);
	CHECKSTAT("aBorderBehavior");
	status = eAttr.addField("Pin", (UCHAR)B::PIN);
	CHECKSTAT("aBorderBehavior");
	status = eAttr.addField("Slide", (UCHAR)B::SLIDE);
	CHECKSTAT("aBorderBehavior");
	status = addAttribute(aBorderBehavior);
	CHECKSTAT("aBorderBehavior");
	status = attributeAffects(aBorderBehavior, outputGeom);
	CHECKSTAT("aBorderBehavior");

	aHardEdgeBehavior = eAttr.create("hardEdgeBehavior", "hb", (UCHAR)B::SLIDE, &status);
	CHECKSTAT("aHardEdgeBehavior");
	eAttr.setKeyable(false);
	eAttr.setChannelBox(true);
	status = eAttr.addField("None", (UCHAR)B::NONE);
	CHECKSTAT("aHardEdgeBehavior");
	status = eAttr.addField("Pin", (UCHAR)B::PIN);
	CHECKSTAT("aHardEdgeBehavior");
	status = eAttr.addField("Slide", (UCHAR)B::SLIDE);
	CHECKSTAT("aHardEdgeBehavior");
	status = addAttribute(aHardEdgeBehavior);
	CHECKSTAT("aHardEdgeBehavior");
	status = attributeAffects(aHardEdgeBehavior, outputGeom);
	CHECKSTAT("aHardEdgeBehavior");

	aGroupEdgeBehavior = eAttr.create("groupEdgeBehavior", "gb", (UCHAR)B::PIN, &status);
	CHECKSTAT("aGroupEdgeBehavior");
	eAttr.setKeyable(false);
	eAttr.setChannelBox(true);
	status = eAttr.addField("None", (UCHAR)B::NONE);
	CHECKSTAT("aGroupEdgeBehavior");
	status = eAttr.addField("Pin", (UCHAR)B::PIN);
	CHECKSTAT("aGroupEdgeBehavior");
	status = eAttr.addField("Slide", (UCHAR)B::SLIDE);
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

	aReprojectDivs = nAttr.create("reprojectDivs", "rpd", MFnNumericData::kInt, 1, &status);
	CHECKSTAT("aReprojectDivs");
	nAttr.setMin(0);
	nAttr.setSoftMax(3);
	nAttr.setChannelBox(true);
	status = addAttribute(aReprojectDivs);
	CHECKSTAT("aReprojectDivs");
	status = attributeAffects(aReprojectDivs, outputGeom);
	CHECKSTAT("aReprojectDivs");

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

	aRecomputeTopo = nAttr.create("recompute", "r", MFnNumericData::kBoolean, false, &status);
	CHECKSTAT("aRecomputeTopo");
	nAttr.setKeyable(false);
	nAttr.setChannelBox(true);
	status = addAttribute(aRecomputeTopo);
	CHECKSTAT("aRecomputeTopo");
	status = attributeAffects(aRecomputeTopo, outputGeom);
	CHECKSTAT("aRecomputeTopo");

	aDeltaMush = nAttr.create("delta", "d", MFnNumericData::kBoolean, false, &status);
	CHECKSTAT("aDeltaMush");
	nAttr.setKeyable(false);
	nAttr.setChannelBox(true);
	status = addAttribute(aDeltaMush);
	CHECKSTAT("aDeltaMush");
	status = attributeAffects(aDeltaMush, outputGeom);
	CHECKSTAT("aDeltaMush");

	aDeltaMult = nAttr.create("deltaMultiplier", "dm", MFnNumericData::kFloat, 1.0f, &status);
	CHECKSTAT("aDeltaMult");
	nAttr.setChannelBox(true);
	status = addAttribute(aDeltaMult);
	CHECKSTAT("aDeltaMult");
	status = attributeAffects(aDeltaMult, outputGeom);
	CHECKSTAT("aDeltaMult");

	aDeltaBase = tAttr.create("deltaBase", "db", MFnData::kMesh);
	tAttr.setArray(true);
	CHECKSTAT("aDeltaBase");
	status = addAttribute(aDeltaBase);
	CHECKSTAT("aDeltaBase");
	status = attributeAffects(aDeltaBase, outputGeom);
	CHECKSTAT("aDeltaBase");


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

		MDataHandle hReprojectDivs = dataBlock.inputValue(aReprojectDivs);
		int reproDivs = hReprojectDivs.asInt();

		MDataHandle hRecompute = dataBlock.inputValue(aRecomputeTopo);
		bool recompute = hRecompute.asBool();

		MDataHandle hTBias = dataBlock.inputValue(aTaubinBias);
		float tBias = hTBias.asFloat();

		MDataHandle hDoDelta = dataBlock.inputValue(aDeltaMush);
		bool doDelta = hDoDelta.asBool();
		MDataHandle hDeltaMult = dataBlock.inputValue(aDeltaMult);
		float deltaMult = hDeltaMult.asFloat();

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

		// Maya automatically copies the input plug to the output plug
		// and then gives you an iterator over that
		// So get the OUTPUT mesh corresponding to the mutiIndex
		MPlug outPlug(thisMObject(), outputGeom);
		outPlug.selectAncestorLogicalIndex(multiIndex, outputGeom);
		MDataHandle hOutput = dataBlock.outputValue(outPlug);
		MObject mesh = hOutput.asMesh();

		// Get the point values
		MFnMesh meshFn(mesh);
		UINT tNumVerts = meshFn.numVertices();
		UINT tNumPolys = meshFn.numPolygons();
		UINT tNumEdges = meshFn.numEdges();

		bool slide = (bb == (UCHAR)B::SLIDE) || (hb == (UCHAR)B::SLIDE) || (gb == (UCHAR)B::SLIDE);

		if (recompute ||
			(bbCheck != bb) || (hbCheck != hb) || (gbCheck != gb) ||
			(tNumVerts != hNumVerts) || (tNumPolys != hNumPolys) || (tNumEdges != hNumEdges)
		) {
			bbCheck = bb; hbCheck = hb; gbCheck = gb;
			hNumVerts = tNumVerts; hNumPolys = tNumPolys; hNumEdges = tNumEdges;
			// Populate the variables with *SPECIALLY ORDERED* data
			// all vertex data is now shuffled by the order vector

			if (relaxer != nullptr) delete relaxer;
			relaxer = MayaRelaxer::Create(mesh, meshFn, vertIter, bb, hb, gb);
		}

		// This can happen if the user is pinning all the points
		// or all the edges are hard (like when you import an obj)
		if (relaxer == nullptr) return status;
		if (relaxer->neighbors.empty()) return status;


		// TODO!! Move these offsets (and the doDelta) into someplace more optimized
		if (doDelta) {

			// get the input mesh corresponding to this output
			MPlug pDeltaBase(thisMObject(), aDeltaBase);
			status = pDeltaBase.selectAncestorLogicalIndex(multiIndex, aDeltaBase);
			CHECKSTAT("deltaBase Select")
			MDataHandle hDeltaBase = dataBlock.inputValue(pDeltaBase, &status);
			if (!status) {
				MGlobal::displayError("Invalid delta base mesh");
				return status;
			}
			MObject deltaBase = hDeltaBase.asMesh();
			if (deltaBase.isNull()) {
				MGlobal::displayError("Null delta base mesh");
				return MStatus::kFailure;
			}

			// Create a copy of the input mesh so we can have maya do the normal/tangent/binormal calculations
			MFnMeshData baseCopyParFn;
			MObject baseCopyPar = baseCopyParFn.create();
			MFnMesh baseCopyFn;
			MObject baseCopy = baseCopyFn.copy(deltaBase, baseCopyPar);

			// get a base and smoothed copy of the original verts
			pointArray_t baseVerts, smoothedBaseVerts;
			baseCopyFn.getPoints(baseVerts);
			bool deltaRepro = false;
			UINT deltaReproDivs = 0;
			smoothedBaseVerts = relaxer->quickRelax(baseCopy, slide, deltaRepro, deltaReproDivs, tBias, iterations);
			baseCopyFn.setPoints(smoothedBaseVerts);
			// TODO: Get the points in n/t/b space and store

			std::vector<MMatrix> mats = relaxer->getVertMatrices(baseCopy, baseCopyFn);
			MItMeshVertex normIt(baseCopy);
			deltas.resize(tNumVerts);
			for (; !normIt.isDone(); normIt.next()) {
				UINT idx = normIt.index();
				deltas[idx] = baseVerts[idx] * mats[idx].inverse();
			}
		}
		else {
			deltas.resize(0);
		}

		// Get the painted weight values
		std::vector<float> weightVals;
		status = getTrueWeights(mesh, dataBlock, multiIndex, weightVals, env);

		// Get the smoothed point positions
		pointArray_t mpa = relaxer->quickRelax(mesh, slide, reproject, reproDivs, tBias, iterations);


		if (doDelta) {
			if (deltas.size() != tNumVerts)
				return MStatus::kFailure;

			for (; !vertIter.isDone(); vertIter.next()) {
				vertIter.setPosition(mpa[vertIter.index()]);
			}
			vertIter.reset();

			std::vector<MMatrix> mats = relaxer->getVertMatrices(mesh, meshFn);

			for (; !vertIter.isDone(); vertIter.next()) {
				int idx = vertIter.index();
				point_t mp = (deltas[idx] * deltaMult) * mats[idx];
				// then do the setPosition on this point
				vertIter.setPosition((weightVals[idx]) * mp + (1.0 - weightVals[idx]) * vertIter.position());
			}
		}
		else {
			for (; !vertIter.isDone(); vertIter.next()) {
				int idx = vertIter.index();
				vertIter.setPosition((weightVals[idx]) * mpa[idx] + (1.0 - weightVals[idx]) * vertIter.position());
			}
		}
	}
	return status;
}
