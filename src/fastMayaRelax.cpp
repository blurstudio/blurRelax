#include <maya/MFnMeshData.h>
#include <maya/MItMeshEdge.h>
#include <maya/MItMeshVertex.h>
#include <maya/MIntArray.h>

#include "fastMayaRelax.h"


MayaRelaxer* MayaRelaxer::Create(
	MObject &mesh, MFnMesh &meshFn, MItGeometry& vertIter,
	short borderBehavior, short hardEdgeBehavior, short groupEdgeBehavior
){ 

	std::vector<std::vector<size_t>> neighbors;
	std::vector<std::vector<UCHAR>> hardEdges;
	std::vector<UCHAR> vertData;

	UINT numVertices = meshFn.numVertices();
	vertData.resize(numVertices);
	hardEdges.resize(numVertices);
	neighbors.resize(numVertices);

	for (; !vertIter.isDone(); vertIter.next()) {
		vertData[vertIter.index()] = (UCHAR)V::IN_GROUP;
	}
	vertIter.reset();

	MItMeshEdge edgeIter(mesh);
	for (; !edgeIter.isDone(); edgeIter.next()) {
		const UINT start = edgeIter.index(0);
		const UINT end = edgeIter.index(1);

		UCHAR edgeData = 0;
		if (edgeIter.onBoundary()) {
			edgeData |= (UCHAR)E::MESH_BORDER;
			vertData[start] |= (UCHAR)V::MESH_BORDER;
			vertData[end] |= (UCHAR)V::MESH_BORDER;
		}
		if (!edgeIter.isSmooth()) edgeData |= (UCHAR)E::HARD;

		neighbors[start].push_back(end);
		neighbors[end].push_back(start);


		if (vertData[start] & (UCHAR)V::IN_GROUP) {
			if (!(vertData[end] & (UCHAR)V::IN_GROUP)) vertData[start] |= (UCHAR)V::GROUP_BORDER;
		}
		else if (vertData[end] & (UCHAR)V::IN_GROUP) vertData[end] |= (UCHAR)V::GROUP_BORDER;

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
				if (!(vertData[polyVerts[i]] & (UCHAR)V::IN_GROUP)){
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

	return new MayaRelaxer(borderBehavior, hardEdgeBehavior, groupEdgeBehavior, neighbors, hardEdges, vertData);
}

void MayaRelaxer::buildOctree(MObject &mesh, bool slide, UINT divisions){
	MStatus status;

	if (divisions == 0){
		octree.create(mesh);
	}
	else {
		MFnMesh meshFn(mesh);
		MFnMeshData smoothMeshParFn;
		MMeshSmoothOptions smoothOpt;
		smoothMeshPar = smoothMeshParFn.create();
		smoothOpt.setDivisions(1);
		smoothOpt.setKeepBorderEdge(slide);
		smoothOpt.setSubdivisionType(MMeshSmoothOptions::kCatmullClark);
		smoothMesh = meshFn.generateSmoothMesh(smoothMeshPar, &smoothOpt);
		octree.create(smoothMesh);
	}
}

void MayaRelaxer::reprojectVerts(FLOAT(*verts)[NUM_COMPS]) const{
	if (!octree.isCreated()) return;
    #pragma omp parallel for if(numVertices>2000)
	for (UINT i = 0; i < numUnpinned; ++i) {
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



MMatrix MayaRelaxer::getMatrixAtPoint(MObject &mesh, MFnMesh &meshFn, MItMeshVertex &vertIt) const {
	MVector norm, tang, binm;
	MIntArray connFaces;


	int idx = vertIt.index();
	vertIt.getNormal(norm);
	vertIt.getConnectedFaces(connFaces);
	meshFn.getFaceVertexTangent(connFaces[0], idx, tang);
	binm = (norm ^ tang).normal();
	tang = (binm ^ norm).normal();
	MPoint vert = vertIt.position();

	/*
	double mm[4][4];
	mm[0][0] = norm[0]; mm[0][1] = tang[0]; mm[0][2] = binm[0]; mm[0][3] = 0.0; 
	mm[1][0] = norm[1]; mm[1][1] = tang[1]; mm[1][2] = binm[1]; mm[1][3] = 0.0;
	mm[2][0] = norm[2]; mm[2][1] = tang[2]; mm[2][2] = binm[2]; mm[2][3] = 0.0;
	mm[3][0] = vert[0]; mm[3][1] = vert[1]; mm[3][2] = vert[2]; mm[3][3] = 1.0;
	*/

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

pointArray_t MayaRelaxer::quickRelax(
	MObject &mesh,
	const bool slide,
	const bool doReproject,
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
		memcpy(&(baseVerts[0][0]), &(verts[0][0]), NUM_COMPS * numVertices * sizeof(FLOAT));
	}

	// Reserve some memory for the next-to-last iteration so we can blend
	FLOAT(*prevVerts)[NUM_COMPS] = new FLOAT[numVertices][NUM_COMPS];

	// If reprojecting, then build the octree
	if (doReproject) buildOctree(mesh, slide, 1u);

	for (size_t r = 0; r < iterI; ++r) {
		// Store the next-to-last iteration to interpolate with
		if ((r == iterI - 1) && (iterT > 0.0))
			memcpy(&(prevVerts[0][0]), &(verts[0][0]), NUM_COMPS * numVertices * sizeof(FLOAT));

		quickLaplacianSmooth(verts);
		if (taubinBias < 1.0)
			quickLaplacianSmooth(verts, taubinBias);
		if (slide)
			edgeProject(baseVerts, verts);
		if (doReproject)
			reprojectVerts(verts);
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

