#include <maya/MFnMeshData.h>
#include <maya/MItMeshEdge.h>
#include <maya/MIntArray.h>

#include "fastMayaRelax.h"


MayaRelaxer MayaRelaxer::Create(
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

	return MayaRelaxer(borderBehavior, hardEdgeBehavior, groupEdgeBehavior, neighbors, hardEdges, vertData);
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

void MayaRelaxer::reprojectVerts(FLOAT(*verts)[4]) const{
	if (!octree.isCreated()) return;
    #pragma omp parallel for if(numVertices>2000)
	for (int i = 0; i < numUnpinned; ++i) {
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

void MayaRelaxer::reorderVerts(MObject &mesh, MFnMesh &meshFn, FLOAT(*reoVerts)[4])const{
	// Build the raw float data buffers
	pointArray_t mpa;
	FLOAT(*reoVerts)[4] = new FLOAT[numVertices][4];
	meshFn.getPoints(mpa);

	for (size_t i = 0; i < numVertices; ++i) {
		reoVerts[i][0] = mpa[(UINT)order[i]].x;
		reoVerts[i][1] = mpa[(UINT)order[i]].y;
		reoVerts[i][2] = mpa[(UINT)order[i]].z;
		reoVerts[i][3] = mpa[(UINT)order[i]].w;
	}
}


void MayaRelaxer::quickRelax(
	MObject &mesh,
	const bool slide,
	const bool doReproject,
	const float taubinBias,
	const FLOAT iterations,
	FLOAT(*verts)[4]
) {
	FLOAT (*baseVerts)[4];
	if (slide) {
		// make a copy of the original verts only if they'll be used for edge reprojection
		baseVerts = new FLOAT[numVertices][4];
		memcpy(&(baseVerts[0][0]), &(verts[0][0]), 4 * numVertices * sizeof(FLOAT));
	}

	FLOAT(*prevVerts)[4];
	prevVerts = new FLOAT[numVertices][4];

	FLOAT iterT, iterFI;
	iterT = modf(iterations, &iterFI);
	UINT iterI = (UINT)iterFI;

	if (iterT > 0.0) {
		iterI += 1;
	}
	
	if (doReproject) buildOctree(mesh, slide, 1u);

	for (size_t r = 0; r < iterI; ++r) {
		// Store the next-to-last iteration to interpolate with
		if ((r == iterI - 1) && (iterT > 0.0))
			memcpy(&(prevVerts[0][0]), &(verts[0][0]), 4 * numVertices * sizeof(FLOAT));

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
		for (size_t i = 0; i < numVertices * 4; ++i) {
			vv[i] = ((vv[i] - pv[i]) * iterT) + pv[i];
		}
	}

	if (slide)
		delete[] baseVerts;
	delete[] prevVerts;
}

