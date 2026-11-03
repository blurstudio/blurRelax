
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
#include <limits>
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
	// creaseIdxs holds only the verts with exactly 2 sliding edges.
	// Those verts have had all their non-sliding neighbors removed, so every
	// stored neighbor is along a hard edge.
	// TODO: Find the hard edge strings and reproject onto those
	for (size_t idx : creaseIdxs) {
		FLOAT *avg = smoothPoints[idx];
		const FLOAT *basePos = basePoints[idx];
		FLOAT keep[3], delta[3], edge[3];

		FLOAT minLen = std::numeric_limits<FLOAT>::max();
		bool found = false;

		// avg - basePos
		for (size_t x = 0; x < 3; ++x)
			delta[x] = avg[x] - basePos[x];

		for (UINT k = neighborOffsets[idx]; k < neighborOffsets[idx + 1]; ++k) {
			const FLOAT *nPos = basePoints[neighborFlat[k]];

			// normalized(basePoints[n] - basePos)
			for (size_t x = 0; x < 3; ++x)
				edge[x] = nPos[x] - basePos[x];
			FLOAT elen2 = edge[0] * edge[0] + edge[1] * edge[1] + edge[2] * edge[2];
			if (elen2 <= 0.0) continue; // degenerate edge

			// dot(delta, edge) / |edge|^2  is the projection length along the unnormalized edge
			FLOAT dd = (delta[0] * edge[0] + delta[1] * edge[1] + delta[2] * edge[2]) / elen2;

			FLOAT dn[3] = { 0.0, 0.0, 0.0 };
			if (dd > 0.0) {
				for (size_t x = 0; x < 3; ++x)
					dn[x] = edge[x] * dd;
			}

			// |delta - dn|^2
			FLOAT len2 = 0.0;
			for (size_t x = 0; x < 3; ++x) {
				FLOAT xx = delta[x] - dn[x];
				len2 += xx * xx;
			}

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
	The verts are sorted by valence (high to low) so all the pinned verts
	(valence 0) are at the end, and we only have to touch the first numUnpinned verts.

	Neighbors are stored in CSR form (neighborOffsets/neighborFlat) so each vert
	can sum its neighbors into registers and write its output exactly once.
	The NUM_COMPS inner loops are fixed-size so they vectorize (4 doubles == one AVX register)

	The update is:
		out = shift * bias * (sum / valence - v) + v
	which is rearranged to avoid the division:
		out = (bias * shift / valence) * sum + (1 - bias * shift) * v
*/
void Relaxer::quickLaplacianSmooth(
	FLOAT verts2d[][NUM_COMPS],
	const FLOAT taubinBias
) const {
	const size_t nzv = numUnpinned;
	if (nzv == 0) return;

	scratch.resize(NUM_COMPS * nzv);
	FLOAT (* __restrict outComp)[NUM_COMPS] = reinterpret_cast<FLOAT(*)[NUM_COMPS]>(scratch.data());
	const FLOAT (* __restrict verts)[NUM_COMPS] = verts2d;
	const UINT * __restrict offsets = neighborOffsets.data();
	const UINT * __restrict flat = neighborFlat.data();
	const FLOAT * __restrict sov = shiftOverValence.data();
	const FLOAT * __restrict shift = shiftVal.data();

	for (size_t i = 0; i < nzv; ++i) {
		FLOAT sum[NUM_COMPS] = {};
		for (UINT k = offsets[i]; k < offsets[i + 1]; ++k) {
			const FLOAT *nv = verts[flat[k]];
			for (size_t j = 0; j < NUM_COMPS; ++j)
				sum[j] += nv[j];
		}
		const FLOAT a = taubinBias * sov[i];
		const FLOAT b = 1.0 - taubinBias * shift[i];
		for (size_t j = 0; j < NUM_COMPS; ++j)
			outComp[i][j] = a * sum[j] + b * verts[i][j];
	}

	std::memcpy(&(verts2d[0][0]), outComp, NUM_COMPS * nzv * sizeof(FLOAT));
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
	std::vector<std::vector<char>> rawSlideEdges; // per-neighbor: does this edge slide
	rawPinPoints.resize(numVertices);
	rawCreaseCount.resize(numVertices);
	rawSlideEdges.resize(numVertices);
	for (size_t start=0; start<rawNeighbors.size(); ++start){
		rawSlideEdges[start].resize(rawNeighbors[start].size(), 0);
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
				else if (slide) { ++rawCreaseCount[start]; rawSlideEdges[start][j] = 1; }
			}
		}
	}

	std::vector<FLOAT> rawShiftVal;
	rawShiftVal.resize(numVertices);
	std::fill(rawShiftVal.begin(), rawShiftVal.end(), 0.5);

	// Apply the "behavior transformations" to the input data
	// if a vert has sliding neighbors, remove all non-sliding neighbors
	// if a vert is pinned remove all of its neighbors
	for (size_t i = 0; i < rawNeighbors.size(); ++i) {
		if ((rawCreaseCount[i] != 0) || rawPinPoints[i]) {
			if (rawCreaseCount[i] != 2)
				rawPinPoints[i] = true;

			std::vector<size_t> newNeigh;
			std::vector<UCHAR> newHard;
			if (!rawPinPoints[i]) {
				for (size_t j = 0; j < rawNeighbors[i].size(); ++j) {
					if (rawSlideEdges[i][j]) {
						newNeigh.push_back(rawNeighbors[i][j]);
						newHard.push_back(rawHardEdges[i][j]);
					}
				}
			}
			rawNeighbors[i] = std::move(newNeigh);
			rawHardEdges[i] = std::move(newHard);
			rawShiftVal[i] = 0.25;
		}
	}

	// Reading the data is done. Now we transform the data for fast processing
	order.resize(numVertices);
	invOrder.resize(numVertices);
	creaseCount.resize(numVertices);
	vertData.resize(numVertices);
	groupIdxs.clear();
	groupIdxs.reserve(numVertices);
	creaseIdxs.clear();
	shiftVal.resize(numVertices);
	valence.resize(numVertices);
	shiftOverValence.resize(numVertices);
	neighborOffsets.clear();
	neighborFlat.clear();
	numUnpinned = 0;

	if (numVertices == 0) return;

	// Sort the verts by their valence count. Compare using > for a sort from high to low
	// Use a stable sort so verts with the same valence keep their original relative order.
	// Mesh vertex order usually has decent spatial locality, and keeping it
	// makes the neighbor lookups in quickLaplacianSmooth much more cache friendly
	// Then build the inverse sort indexes
	std::iota(order.begin(), order.end(), 0u);
	std::stable_sort(order.begin(), order.end(), [&rawNeighbors](size_t a, size_t b) {return rawNeighbors[a].size() > rawNeighbors[b].size(); });
	for (size_t i = 0; i < order.size(); ++i) { invOrder[order[i]] = i; }

	size_t totalNeighbors = 0;
	for (const auto &neigh : rawNeighbors) totalNeighbors += neigh.size();
	neighborFlat.reserve(totalNeighbors);
	neighborOffsets.reserve(numVertices + 1);
	neighborOffsets.push_back(0);

	// Build the CSR neighbor lists
	for (size_t i = 0; i < numVertices; ++i) {
		const auto &neigh = rawNeighbors[order[i]];
		size_t vale = neigh.size();
		if (vale > 0) {
			numUnpinned = (UINT)(i + 1);
			for (size_t n = 0; n < vale; ++n) {
				neighborFlat.push_back((UINT)invOrder[neigh[n]]);
			}
			neighborOffsets.push_back((UINT)neighborFlat.size());
		}
		// only needed for slide and reproject
		// and I only need to know if it's ==0 or ==2
		creaseCount[i] = rawCreaseCount[order[i]];
		vertData[i] = rawVertData[order[i]];

		if (vertData[i] & (UCHAR)V::IN_GROUP)
			groupIdxs.push_back((size_t)i);
		if (creaseCount[i] == 2 && vale > 0)
			creaseIdxs.push_back((size_t)i);

		valence[i] = FLOAT(vale);
		shiftVal[i] = rawShiftVal[order[i]];
		shiftOverValence[i] = (vale > 0) ? shiftVal[i] / valence[i] : 0.0;
	}
}



