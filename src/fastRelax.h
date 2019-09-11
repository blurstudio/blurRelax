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

typedef double FLOAT;
typedef unsigned int UINT;
#define NUM_COMPS 4

// Yes, these should be enums. I'll get to it
#define V_IN_GROUP 1
#define V_MESH_BORDER 2
#define V_GROUP_BORDER 4

#define E_HARD 1
#define E_MESH_BORDER 2
#define E_GROUP_BORDER 4

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


void edgeProject(
	const FLOAT basePoints[][NUM_COMPS],
	const std::vector<size_t> &group,
	const std::vector<size_t> &invOrder,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<UINT> &creaseCount,
	FLOAT smoothPoints[][NUM_COMPS]
);

void quickLaplacianSmooth(
	FLOAT verts2d[][NUM_COMPS],
	const size_t numVerts,
	const std::vector<std::vector<size_t>> &neighbors,
	const std::vector<FLOAT> &valence,
	const std::vector<FLOAT> &shiftVal,
	const FLOAT taubinBias=1.0
);

void fillQuickTopoVars(
	// Behaviors
	short borderBehavior, // BB_NONE/BB_PIN/BB_SLIDE
	short hardEdgeBehavior, // HB_NONE/HB_PIN/HB_SLIDE
	short groupEdgeBehavior, // GB_NONE/GB_PIN/GB_SLIDE

	// Inputs
	std::vector<std::vector<size_t>> rawNeighbors, // A vector of neighbor indices per vertex. Copied
	std::vector<std::vector<char>> rawHardEdges, // Bitwise per-neighbor data: edge is hard, edge along boundary. Copied
	const std::vector<char> &rawVertData, // Bitwise per-vert data: Group membership, geo boundary, group boundary,

	// Outputs
	std::vector<std::vector<size_t>> &neighbors,
	std::vector<UINT> &creaseCount,
	std::vector<FLOAT> &shiftVal,
	std::vector<FLOAT> &valence,
	std::vector<size_t> &order,
	std::vector<size_t> &invOrder
);
