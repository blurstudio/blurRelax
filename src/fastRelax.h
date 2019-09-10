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
#pragma once
#include <vector>
#include <algorithm>
#include <numeric>

#define V_IN_GROUP 1
#define V_MESH_BORDER 2
#define V_GROUP_BORDER 4

#define E_HARD 1
#define E_MESH_BORDER 2
#define E_GROUP_BORDER 4

#define NUM_COMPS 4


void edgeProject(
	const float_t basePoints[][NUM_COMPS],
	const std::vector<size_t> &group,
	const std::vector<size_t> &invOrder,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<UINT> &creaseCount,
	float_t smoothPoints[][NUM_COMPS]
);

void quickLaplacianSmooth(
	float_t verts2d[][NUM_COMPS],
	const size_t numVerts,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<float_t> &valence,
	const std::vector<float_t> &shiftVal,
	const float_t taubinBias=1.0
);

void loadMayaTopologyData(
	MObject &mesh,
	MItGeometry& vertIter,
	std::vector<std::vector<UINT>> &neighbors, // A vector of neighbor indices per vertex
	std::vector<std::vector<char>> &hardEdges, // Bitwise per-neighbor data: edge is hard, edge along boundary
	std::vector<char> &vertData // Bitwise per-vert data: Group membership, geo boundary, group boundary,
);

void fillQuickTopoVars(
	// Behaviors
	short borderBehavior, // BB_NONE/BB_PIN/BB_SLIDE
	short hardEdgeBehavior, // HB_NONE/HB_PIN/HB_SLIDE
	short groupEdgeBehavior, // GB_NONE/GB_PIN/GB_SLIDE

	// Inputs
	std::vector<std::vector<UINT>> rawNeighbors, // A vector of neighbor indices per vertex. Copied
	std::vector<std::vector<char>> rawHardEdges, // Bitwise per-neighbor data: edge is hard, edge along boundary. Copied
	const std::vector<char> &rawVertData // Bitwise per-vert data: Group membership, geo boundary, group boundary,

	// Outputs
	std::vector<std::vector<size_t>> &neighbors,
	std::vector<UINT> &creaseCount,
	std::vector<float_t> &shiftVal,
	std::vector<float_t> &valence,
	std::vector<size_t> &order,
	std::vector<size_t> &invOrder
);

