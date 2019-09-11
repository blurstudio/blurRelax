
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
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include "fastRelax.h"
#include "blurRelaxNode.h"

/*
	This code takes the edges that have been defined as hard
	and reprojects them onto the original neighboring hard edges

	Perhaps in the future, I can project onto non-branching
	strings of edges. But for now, it's just direct neighbors
*/


void edgeProject(
	const float_t basePoints[][4],
	const std::vector<size_t> &group,
	const std::vector<size_t> &invOrder,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<UINT> &creaseCount,
	float_t smoothPoints[][4]
) {
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


/*
	All the crazy hoops I've jumped through are to make auto-vectorization work
	in this function.

	neighbors and neighborOffsets work together to be a kind of pseudo-transpose
	of a std::vector of neighbors per vert with that vector sorted so the verts
	with the most neighbors were at the top.
	neighborOffsets contains the index offsets of vertices with at least [index] of neighbors
*/
void quickLaplacianSmooth(
	float_t verts2d[][4],
	const size_t numVerts,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<float_t> &valence,
	const std::vector<float_t> &shiftVal,
	const float_t taubinBias=1.0
) {
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


/*
	I'm trying to pre-process everything I can right here so I don't have to
	have any 'if's in quickLaplacianSmooth() so auto-vectorization works
*/
void fillQuickTopoVars(
	// Behaviors
	short borderBehavior, // BB_NONE/BB_PIN/BB_SLIDE
	short hardEdgeBehavior, // HB_NONE/HB_PIN/HB_SLIDE
	short groupEdgeBehavior, // GB_NONE/GB_PIN/GB_SLIDE

	// Inputs
	std::vector<std::vector<UINT>> rawNeighbors, // A vector of neighbor indices per vertex. Copied
	std::vector<std::vector<char>> rawHardEdges, // Bitwise per-neighbor data: edge is hard, edge along boundary. Copied
	const std::vector<char> &rawVertData, // Bitwise per-vert data: Group membership, geo boundary, group boundary,

	// Outputs
	std::vector<std::vector<size_t>> &neighbors,
	std::vector<UINT> &creaseCount,
	std::vector<float_t> &shiftVal,
	std::vector<float_t> &valence,
	std::vector<size_t> &order,
	std::vector<size_t> &invOrder
){
	size_t numVertices = rawVertData.size();

	// Read the input data and define the per-point behavior
	std::vector<char> rawPinPoints;
	std::vector<UINT> rawCreaseCount;
	rawPinPoints.resize(numVertices);
	rawCreaseCount.resize(numVertices);
	for (size_t start=0; start<rawNeighbors.size(); ++start){
		const char startData = rawVertData[start];
		if (!(startData & V_IN_GROUP)) continue; // short circuit if the start isn't in the group
		for (size_t j=0; j<rawNeighbors[start].size(); ++j){
			const UINT end = rawNeighbors[start][j];
			const char endData = rawVertData[end];
			if (endData & V_IN_GROUP){ // start and end are both in group
				const char edgeData = rawHardEdges[start][j];

				// Figure out if this edge is pinned
				const bool pin =
					((borderBehavior    == BB_PIN) && (edgeData & E_MESH_BORDER)) ||
					((hardEdgeBehavior  == HB_PIN) && (edgeData & E_HARD)         && !(edgeData & E_MESH_BORDER)) ||
					((groupEdgeBehavior == GB_PIN) && (edgeData & E_GROUP_BORDER) && !(edgeData & E_MESH_BORDER));

				// Figure out if this edge slides
				const bool slide =
					((borderBehavior    == BB_SLIDE) && (edgeData & E_MESH_BORDER)) ||
					((hardEdgeBehavior  == HB_SLIDE) && (edgeData & E_HARD)         && !(edgeData & E_MESH_BORDER)) ||
					((groupEdgeBehavior == GB_SLIDE) && (edgeData & E_GROUP_BORDER) && !(edgeData & E_MESH_BORDER));

				if (pin) { rawPinPoints[start] = true; }
				else if (slide) { ++rawCreaseCount[start]; }
			}
		}
	}

	std::vector<float_t> rawShiftVal;
	rawShiftVal.resize(numVertices);
	std::fill(rawShiftVal.begin(), rawShiftVal.end(), 0.5);

	// Apply the "behavior transformations" to the input data
	// if a vert has hard neighbors, remove all non-hard neighbors
	// if a vert is pinned remove all of its neighbors
	for (size_t i = 0; i < rawNeighbors.size(); ++i) {
		bool pinThisPoint = rawPinPoints[i];
		if ((rawCreaseCount[i] != 0) || pinThisPoint) {
			rawShiftVal[i] = 0.25;
			if (rawCreaseCount[i] != 2) pinThisPoint = true;
			if (!pinThisPoint) {
				std::vector<UINT> newNeigh;
				std::vector<char> newHard;
				for (size_t j = 0; j < rawNeighbors[i].size(); ++j) {
					if (rawHardEdges[i][j] & E_HARD) {
						newNeigh.push_back(rawNeighbors[i][j]);
						newHard.push_back(rawHardEdges[i][j]);
					}
				}
				rawNeighbors[i] = newNeigh;
				rawHardEdges[i] = newHard;
			}
		}
	}

	// Reading the data is done. Now we transform the data for fast processing
	size_t maxValence = rawNeighbors[order[0]].size();
	neighbors.resize(maxValence);
	order.resize(numVertices);
	invOrder.resize(numVertices);
	creaseCount.resize(numVertices);
	shiftVal.resize(numVertices * NUM_COMPS);
	valence.resize(numVertices * NUM_COMPS);

	// Sort the verts by their valence count. Compare using > for a sort from high to low
	// Then build the inverse sort indexes
	std::iota(order.begin(), order.end(), 0u);
	std::sort(order.begin(), order.end(), [&rawNeighbors](size_t a, size_t b) {return rawNeighbors[a].size() > rawNeighbors[b].size(); });
	for (size_t i = 0; i < order.size(); ++i) { invOrder[order[i]] = i; }

	// Build the "transposed" neighbor and hard edge values
	for (size_t i = 0; i < numVertices; ++i) {
		const auto &neigh = rawNeighbors[order[i]];
		const auto &hards = rawHardEdges[order[i]];
		size_t vale = neigh.size();
		for (size_t n = 0; n < vale; ++n) {
			neighbors[n].push_back(invOrder[neigh[n]]);
		}
		// only needed for slide and reproject
		// and I only need to know if it's ==0 or ==2
		creaseCount[i] = rawCreaseCount[order[i]];

		// Vectorizing flattens the vert list, so I need this data per vert, per component
		// Maya uses xyzw points, so I need 4. In other cases I'll need 3
		for (size_t xx = 0; xx < NUM_COMPS; ++xx) {
			valence[NUM_COMPS * i + xx] = float_t(vale);
			shiftVal[NUM_COMPS * i + xx] = rawShiftVal[order[i]];
		}
	}
}


/*
   Load the minimal topology data from Maya
*/
void loadMayaTopologyData(
	// Inputs
	MObject &mesh,
	MItGeometry& vertIter,

	//outputs
	std::vector<std::vector<UINT>> &neighbors, // A vector of neighbor indices per vertex
	std::vector<std::vector<char>> &hardEdges, // Bitwise per-neighbor data: edge is hard, edge along boundary
	std::vector<char> &vertData // Bitwise per-vert data: Group membership, geo boundary, group boundary,
){
	MFnMesh meshFn(mesh);
	UINT numVertices = meshFn.numVertices();
	vertData.resize(numVertices);
	hardEdges.resize(numVertices);
	neighbors.resize(numVertices);

	for (; !vertIter.isDone(); vertIter.next()) {
		vertData[vertIter.index()] = V_IN_GROUP;
	}
	vertIter.reset();

	MItMeshEdge edgeIter(mesh);
	for (; !edgeIter.isDone(); edgeIter.next()) {
		const UINT start = edgeIter.index(0);
		const UINT end = edgeIter.index(1);

		char edgeData;
		if (edgeIter.onBoundary()) {
			edgeData |= E_MESH_BORDER;
			vertData[start] |= V_MESH_BORDER;
			vertData[end] |= V_MESH_BORDER;
		}
		if (!edgeIter.isSmooth()) edgeData |= E_HARD;

		neighbors[start].push_back(end);
		neighbors[end].push_back(start);


		if (vertData[start] & V_IN_GROUP) {
			if (!(vertData[end] & V_IN_GROUP)) vertData[start] |= V_GROUP_BORDER;
		}
		else if (vertData[end] & V_IN_GROUP) vertData[end] |= V_GROUP_BORDER;

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
				if (!(vertData[polyVerts[i]] & V_IN_GROUP)){
					internalEdge = false;
					edgeData |= E_GROUP_BORDER;
					break;
				}
			}
			if (!internalEdge) break;
		}

		hardEdges[start].push_back(edgeData);
		hardEdges[end].push_back(edgeData);

	}
}


