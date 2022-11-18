
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
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "fastRelax.h"

/*
	This code takes the edges that have been defined as hard
	and reprojects them onto the original neighboring hard edges

	Perhaps in the future, I can project onto non-branching
	strings of edges. But for now, it's just direct neighbors
*/

void Relaxer::edgeProject(
	const FLOAT basePoints[][NUM_COMPS],
	FLOAT smoothPoints[][NUM_COMPS]
) const {
	std::vector<size_t> neigh;
	for (size_t gidx = 0; gidx < groupIdxs.size(); ++gidx) {
		// If we have "hard edges" we have already removed
		// any soft edges, so if the crease count for this vertex is zero
		// we can just completely skip it. And if it's >0, then we can
		// assume that all stored neighbors are hard.
		size_t idx = invOrder[groupIdxs[gidx]];

		FLOAT *avg = smoothPoints[idx];
		const FLOAT *basePos = basePoints[idx];
		if (creaseCount[idx] != 2) { continue; }
		FLOAT keep[3], delta[3], edge[3];

		neigh.clear();
		for (size_t i = 0; i < neighbors.size(); ++i) {
			if (idx < neighbors[i].size()) {
				neigh.push_back(neighbors[i][idx]);
			}
			else break;
		}

		// TODO: Find the hard edge strings and reproject onto those
		FLOAT minLen = std::numeric_limits<FLOAT>::max();
		bool found = false;

		// avg - basePos
		for (size_t x = 0; x < 3; ++x)
			delta[x] = avg[x] - basePos[x];

		for (size_t i = 0; i < neigh.size(); ++i) {
			size_t n = neigh[i];

			// normalized(prevPoints[n] - basePos)
			for (size_t x = 0; x < 3; ++x)
				edge[x] = basePoints[n][x] - basePos[x];
			FLOAT elen = std::sqrt(edge[0] * edge[0] + edge[1] * edge[1] + edge[2] * edge[2]);
			for (size_t x = 0; x < 3; ++x)
				edge[x] /= elen;

			// dot(delta, edge)
			FLOAT dd = delta[0] * edge[0] + delta[1] * edge[1] + delta[2] * edge[2];

			FLOAT dn[3] = { 0.0f, 0.0f, 0.0f };
			if (dd > 0.0) {
				for (size_t x = 0; x < 3; ++x)
					dn[x] = edge[x] * dd;
			}
			FLOAT xx[3];

			// delta - dn
			for (size_t x = 0; x < 3; ++x)
				xx[x] = delta[x] - dn[x];

			// dot(xx, xx)
			FLOAT len2 = xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2];

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

void Relaxer::quickLaplacianSmooth(
	FLOAT verts2d[][NUM_COMPS],
	const FLOAT taubinBias
) const {
	// First, get verts as a single pointer to the contiguous memory stored in (verts2d*)[4]
	FLOAT* verts = &(verts2d[0][0]);

	// number of nonzero valence
	size_t nzv = neighbors[0].size();

	// number of nonzero components
	size_t nzc = NUM_COMPS * nzv;

	// The __restrict keyword tells the compiler that *outComp
	// is not pointed to by any other pointer in this scope
	// This allows for auto-vectorization
	FLOAT * __restrict outComp = new FLOAT[nzc];
	memset(outComp, 0, nzc*sizeof(FLOAT));

	for (size_t ncIdx = 0; ncIdx < neighbors.size(); ++ncIdx) {
		const auto &nCol = neighbors[ncIdx];
		size_t nColCount = nCol.size();
		for (size_t i = 0; i < nColCount; ++i) {
			size_t nci = NUM_COMPS * nCol[i];
			for (size_t j = 0; j < NUM_COMPS; ++j) {
				outComp[NUM_COMPS*i+j] = outComp[NUM_COMPS*i+j] + verts[nci+j];
			}
		}
	}

	for (size_t i = 0; i < nzc; ++i) {
		outComp[i] = shiftVal[i] * taubinBias * ((outComp[i] / valence[i]) - verts[i]) + verts[i];
	}

	std::memcpy(verts, outComp, nzc*sizeof(FLOAT));
	delete outComp;
}


/*
	I'm trying to pre-process everything I can right here so I don't have to
	have any 'if's in quickLaplacianSmooth() so auto-vectorization works

	borderBehavior, // BB_NONE/BB_PIN/BB_SLIDE
	hardEdgeBehavior, // HB_NONE/HB_PIN/HB_SLIDE
	groupEdgeBehavior, // GB_NONE/GB_PIN/GB_SLIDE
	std::vector<std::vector<size_t>> rawNeighbors, // A vector of neighbor indices per vertex
	std::vector<std::vector<UCHAR>> rawHardEdges, // Bitwise per-neighbor data: edge is hard, edge along boundary
	const std::vector<UCHAR> &rawVertData // Bitwise per-vert data: Group membership, geo boundary, group boundary,

*/
Relaxer::Relaxer(
	short borderBehavior,
	short hardEdgeBehavior,
	short groupEdgeBehavior,
	std::vector<std::vector<size_t>> rawNeighbors, // Not a reference
	std::vector<std::vector<UCHAR>> rawHardEdges, // Not a reference
	const std::vector<UCHAR> &rawVertData
) {

	// Read the input data and define the per-point behavior
	numVertices = (UINT)rawNeighbors.size();
	std::vector<UCHAR> rawPinPoints;
	std::vector<UINT> rawCreaseCount;
	rawPinPoints.resize(numVertices);
	rawCreaseCount.resize(numVertices);
	for (size_t start=0; start<rawNeighbors.size(); ++start){
		const UCHAR startData = rawVertData[start];
		if (!(startData & (UCHAR)V::IN_GROUP)) continue; // short circuit if the start isn't in the group
		for (size_t j=0; j<rawNeighbors[start].size(); ++j){
			const UINT end = (UINT)rawNeighbors[start][j];
			const UCHAR endData = rawVertData[end];
			if (endData & (UCHAR)V::IN_GROUP){ // start and end are both in group
				const UCHAR edgeData = rawHardEdges[start][j];

				// Figure out if this edge is pinned
				const bool pin =
					((borderBehavior    == (UCHAR)B::PIN) && (edgeData & (UCHAR)E::MESH_BORDER)) ||
					((hardEdgeBehavior  == (UCHAR)B::PIN) && (edgeData & (UCHAR)E::HARD)         && !(edgeData & (UCHAR)E::MESH_BORDER)) ||
					((groupEdgeBehavior == (UCHAR)B::PIN) && (edgeData & (UCHAR)E::GROUP_BORDER) && !(edgeData & (UCHAR)E::MESH_BORDER));

				// Figure out if this edge slides
				const bool slide =
					((borderBehavior    == (UCHAR)B::SLIDE) && (edgeData & (UCHAR)E::MESH_BORDER)) ||
					((hardEdgeBehavior  == (UCHAR)B::SLIDE) && (edgeData & (UCHAR)E::HARD)         && !(edgeData & (UCHAR)E::MESH_BORDER)) ||
					((groupEdgeBehavior == (UCHAR)B::SLIDE) && (edgeData & (UCHAR)E::GROUP_BORDER) && !(edgeData & (UCHAR)E::MESH_BORDER));

				if (pin) { rawPinPoints[start] = true; }
				else if (slide) { ++rawCreaseCount[start]; }
			}
		}
	}

	std::vector<FLOAT> rawShiftVal;
	rawShiftVal.resize(numVertices);
	std::fill(rawShiftVal.begin(), rawShiftVal.end(), 0.5);

	// Apply the "behavior transformations" to the input data
	// if a vert has hard neighbors, remove all non-hard neighbors
	// if a vert is pinned remove all of its neighbors
	for (size_t i = 0; i < rawNeighbors.size(); ++i) {
		if ((rawCreaseCount[i] != 0) || rawPinPoints[i]) {
			if (rawCreaseCount[i] != 2)
				rawPinPoints[i] = true;

			std::vector<size_t> newNeigh;
			std::vector<UCHAR> newHard;
			if (!rawPinPoints[i]) {
				for (size_t j = 0; j < rawNeighbors[i].size(); ++j) {
					if (rawHardEdges[i][j] & (UCHAR)E::HARD) {
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

	// Reading the data is done. Now we transform the data for fast processing
	order.resize(numVertices);
	invOrder.resize(numVertices);
	creaseCount.resize(numVertices);
	vertData.resize(numVertices);
	groupIdxs.resize(numVertices);
	shiftVal.resize(numVertices * NUM_COMPS);
	valence.resize(numVertices * NUM_COMPS);

	// Sort the verts by their valence count. Compare using > for a sort from high to low
	// Then build the inverse sort indexes
	std::iota(order.begin(), order.end(), 0u);
	std::sort(order.begin(), order.end(), [&rawNeighbors](size_t a, size_t b) {return rawNeighbors[a].size() > rawNeighbors[b].size(); });
	for (size_t i = 0; i < order.size(); ++i) { invOrder[order[i]] = i; }

	size_t maxValence = rawNeighbors[order[0]].size();
	neighbors.resize(maxValence);

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
		vertData[i] = rawVertData[order[i]];

		if (vertData[i] & (UCHAR)V::IN_GROUP)
			groupIdxs.push_back((size_t)i);

		// Vectorizing flattens the vert list, so I need this data per vert, per component
		// Maya uses xyzw points, so I need 4. In other cases I'll need 3
		for (size_t xx = 0; xx < NUM_COMPS; ++xx) {
			valence[NUM_COMPS * i + xx] = FLOAT(vale);
			shiftVal[NUM_COMPS * i + xx] = rawShiftVal[order[i]];
		}
	}
}


