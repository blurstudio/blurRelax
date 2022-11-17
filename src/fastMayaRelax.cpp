#include <cstring>
#include <maya/MFnMeshData.h>
#include <maya/MItMeshEdge.h>
#include <maya/MItMeshVertex.h>
#include <maya/MIntArray.h>
#include <maya/MItMeshPolygon.h>

#include "fastMayaRelax.h"


MayaRelaxer* MayaRelaxer::Create(
	MObject &mesh, MFnMesh &meshFn, MItGeometry& vertIter,
	short borderBehavior, short hardEdgeBehavior, short groupEdgeBehavior
){ 

	std::vector<std::vector<size_t>> neighbors;
	std::vector<std::vector<UCHAR>> hardEdges;
	std::vector<UCHAR> rawVertData;

	UINT numVertices = meshFn.numVertices();
	rawVertData.resize(numVertices);
	hardEdges.resize(numVertices);
	neighbors.resize(numVertices);

	for (; !vertIter.isDone(); vertIter.next()) {
		rawVertData[vertIter.index()] = (UCHAR)V::IN_GROUP;
	}
	vertIter.reset();

	MItMeshEdge edgeIter(mesh);
	for (; !edgeIter.isDone(); edgeIter.next()) {
		const UINT start = edgeIter.index(0);
		const UINT end = edgeIter.index(1);

		UCHAR edgeData = 0;
		if (edgeIter.onBoundary()) {
			edgeData |= (UCHAR)E::MESH_BORDER;
			rawVertData[start] |= (UCHAR)V::MESH_BORDER;
			rawVertData[end] |= (UCHAR)V::MESH_BORDER;
		}
		if (!edgeIter.isSmooth()) edgeData |= (UCHAR)E::HARD;

		neighbors[start].push_back(end);
		neighbors[end].push_back(start);


		if (rawVertData[start] & (UCHAR)V::IN_GROUP) {
			if (!(rawVertData[end] & (UCHAR)V::IN_GROUP)) rawVertData[start] |= (UCHAR)V::GROUP_BORDER;
		}
		else if (rawVertData[end] & (UCHAR)V::IN_GROUP) rawVertData[end] |= (UCHAR)V::GROUP_BORDER;

		// an edge is a group border edge iff 
		// one of the two faces bordering this edge
		// has a vertex that is not in the group
		bool internalEdge = true;
		MIntArray connFaces;
		edgeIter.getConnectedFaces(connFaces);
		for (UINT i=0; i<connFaces.length(); ++i){
			MIntArray polyVerts;
			meshFn.getPolygonVertices(connFaces[i], polyVerts);
			for (UINT j=0; j<polyVerts.length(); ++j){
				if (!(rawVertData[polyVerts[i]] & (UCHAR)V::IN_GROUP)){
					internalEdge = false;
					edgeData |= (UCHAR)E::GROUP_BORDER;
					break;
				}
			}
			if (!internalEdge) break;
		}

		hardEdges[start].push_back(edgeData);
		hardEdges[end].push_back(edgeData);
	}

	return new MayaRelaxer(borderBehavior, hardEdgeBehavior, groupEdgeBehavior, neighbors, hardEdges, rawVertData);
}

std::vector<MMatrix> MayaRelaxer::getVertMatrices(MObject &mesh, MFnMesh &meshFn) const {
	MStatus stat;
	UINT numVerts = meshFn.numVertices();

	MFloatVectorArray norms, faceTangs, tangs;
	MFloatPointArray points;
	meshFn.getVertexNormals(false, norms);
	meshFn.getTangents(faceTangs);
	meshFn.getPoints(points);

	tangs.setLength(numVerts);
	
	MIntArray counts, connects;
	meshFn.getVertices(counts, connects);

	int offset = 0;
	std::vector<char> dones(numVerts);
	for (UINT faceIdx = 0; faceIdx < counts.length(); ++faceIdx) {
		int count = counts[faceIdx];

		for (int i = 0; i < count; ++i) {
			int vIdx = connects[offset + i];
			if (!dones[vIdx]) {
				dones[vIdx] = 1;
				tangs[vIdx] = faceTangs[offset + i];
			}
		}
		offset += count;
	}

	std::vector<MMatrix> out(numVerts);
	for (UINT i = 0; i < numVerts; ++i) {
		MFloatPoint &vert = points[i];
		MFloatVector &norm = norms[i];
		MFloatVector &tang = tangs[i];

		MFloatVector binm = (norm ^ tang).normal();
		tang = (binm ^ norm).normal();

		double mm[4][4] = {
			{norm[0], tang[0], binm[0], 0.0},
			{norm[1], tang[1], binm[1], 0.0},
			{norm[2], tang[2], binm[2], 0.0},
			{vert[0], vert[1], vert[2], 1.0}
		};
		out[i] = MMatrix(mm);
	}

	return out;
}



MMatrix MayaRelaxer::getMatrixAtPoint(MObject &mesh, MFnMesh &meshFn, MItMeshVertex &vertIt) const {
	MVector norm, tang, binm;

	int idx = vertIt.index();
	vertIt.getNormal(norm);
	meshFn.getFaceVertexTangent(connFaces[idx], idx, tang);
	binm = (norm ^ tang).normal();
	tang = (binm ^ norm).normal();
	MPoint vert = vertIt.position();

	double mm[4][4] = {
		{norm[0], tang[0], binm[0], 0.0},
		{norm[1], tang[1], binm[1], 0.0},
		{norm[2], tang[2], binm[2], 0.0},
		{vert[0], vert[1], vert[2], 1.0}
	};

	MMatrix mat(mm);
	return mat;
}

void MayaRelaxer::reorderVerts(pointArray_t mpa, FLOAT(*reoVerts)[NUM_COMPS])const{
	// Build the raw float data buffers
	for (size_t i = 0; i < numVertices; ++i) {
		for (size_t j = 0; j < NUM_COMPS; ++j) {
			reoVerts[i][j] = mpa[(UINT)order[i]][(UINT)j];
		}
	}
}

void MayaRelaxer::reorderVerts(MObject &mesh, MFnMesh &meshFn, FLOAT(*reoVerts)[NUM_COMPS]) const{
	// Build the raw float data buffers
	pointArray_t mpa;
	meshFn.getPoints(mpa);
	reorderVerts(mpa, reoVerts);
}

void MayaRelaxer::revertVerts(FLOAT(*reoVerts)[NUM_COMPS], FLOAT(*reverted)[NUM_COMPS]) const{
	for (size_t i = 0; i < numVertices; ++i) {
		for (size_t j = 0; j < NUM_COMPS; ++j) {
			reverted[i][j] = reoVerts[(UINT)invOrder[i]][j];
		}
	}
}

pointArray_t MayaRelaxer::revertVerts(FLOAT(*reoVerts)[NUM_COMPS]) const{
	FLOAT(*rawRevert)[NUM_COMPS] = new FLOAT[numVertices][NUM_COMPS];
	revertVerts(reoVerts, rawRevert);
	pointArray_t reverted(rawRevert, numVertices);
	delete[] rawRevert;
	return reverted;
}




void MayaRelaxer::buildOctree(MObject &mesh, bool slide, UINT divisions){
}

void MayaRelaxer::reprojectVerts(FLOAT(*verts)[NUM_COMPS], MMeshIntersector &octree) const{
	if (!octree.isCreated()) return;
	UINT nzv = (UINT)neighbors[0].size();
    #pragma omp parallel for if(numVertices>2000)
	for (UINT i = 0; i < nzv; ++i) {
		if ((creaseCount[i] == 0) && (vertData[i] && (UCHAR)V::IN_GROUP)) {
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





pointArray_t MayaRelaxer::quickRelax(
	MObject &mesh,
	const bool slide,
	const bool doReproject,
	const UINT reprojectDivisons,
	const float taubinBias,
	const FLOAT iterations
) {

	// Split the float iterations into an int and float
	FLOAT iterT, iterFI;
	iterT = modf(iterations, &iterFI);
	UINT iterI = (UINT)iterFI;
	if (iterT > 0.0) { iterI += 1; }

	// Get the mesh verts
	MFnMesh meshFn(mesh);
	FLOAT(*verts)[NUM_COMPS] = new FLOAT[numVertices][NUM_COMPS];
	reorderVerts(mesh, meshFn, verts);

	// make a copy of the original verts only if they'll be used for edge reprojection
	FLOAT(*baseVerts)[NUM_COMPS];
	if (slide) {
		baseVerts = new FLOAT[numVertices][NUM_COMPS];
		std::memcpy(&(baseVerts[0][0]), &(verts[0][0]), NUM_COMPS * numVertices * sizeof(FLOAT));
	}

	// Reserve some memory for the next-to-last iteration so we can blend
	FLOAT(*prevVerts)[NUM_COMPS] = new FLOAT[numVertices][NUM_COMPS];

	// If reprojecting, then build the octree
	MMeshIntersector octree;
	MObject smoothMeshPar, smoothMesh;
	if (doReproject) {
		MStatus status;
		int divisions = 0;

		if (reprojectDivisons == 0){
			octree.create(mesh);
		}
		else {
			MFnMesh meshFn(mesh);
			MFnMeshData smoothMeshParFn;
			MMeshSmoothOptions smoothOpt;
			smoothMeshPar = smoothMeshParFn.create();
			smoothOpt.setDivisions(reprojectDivisons);
			smoothOpt.setKeepBorderEdge(slide);
			smoothOpt.setSubdivisionType(MMeshSmoothOptions::kCatmullClark);
			smoothMesh = meshFn.generateSmoothMesh(smoothMeshPar, &smoothOpt);
			octree.create(smoothMesh);
		}
	}

	for (size_t r = 0; r < iterI; ++r) {
		// Store the next-to-last iteration to interpolate with
		if ((r == iterI - 1) && (iterT > 0.0))
			std::memcpy(&(prevVerts[0][0]), &(verts[0][0]), NUM_COMPS * numVertices * sizeof(FLOAT));

		quickLaplacianSmooth(verts);
		if (taubinBias < 1.0)
			quickLaplacianSmooth(verts, taubinBias);
		if (slide)
			edgeProject(baseVerts, verts);
		if (doReproject)
			reprojectVerts(verts, octree);
	}

	// Interpolate between prevVerts and verts based on iterT
	if (iterT > 0.0) {
		// This should vectorize
		FLOAT * vv = &verts[0][0];
		FLOAT * pv = &prevVerts[0][0];
		for (size_t i = 0; i < numVertices * NUM_COMPS; ++i) {
			vv[i] = ((vv[i] - pv[i]) * iterT) + pv[i];
		}
	}

	pointArray_t ret = revertVerts(verts);

	// clean up after myself
	if (slide)
		delete[] baseVerts;
	delete[] prevVerts;
	delete[] verts;

	return ret;
}

