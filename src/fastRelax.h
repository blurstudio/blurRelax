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
typedef unsigned char UCHAR;
constexpr auto NUM_COMPS = 4;


enum class V : UCHAR { IN_GROUP = 1, MESH_BORDER = 2, GROUP_BORDER = 4 };
enum class E : UCHAR { HARD = 1, MESH_BORDER = 2, GROUP_BORDER = 4 };
enum class B : UCHAR { NONE, PIN, SLIDE };
enum class S : UCHAR { TAUBIN, LAPLACIAN };


class Relaxer {
public:
	std::vector<size_t> order;
	std::vector<size_t> invOrder;
	std::vector<size_t> groupIdxs;
	std::vector<std::vector<size_t>> neighbors; // A vector of neighbor indices per vertex
	std::vector<std::vector<UCHAR>> hardEdges; // Bitwise per-neighbor data: edge is hard, edge along boundary
	std::vector<UCHAR> vertData; // Bitwise per-vert data: Group membership, geo boundary, group boundary,
	std::vector<FLOAT> shiftVal; // normally 0.5; but it's 0.25 if on a hard edge
	std::vector<FLOAT> valence; // as float for vectorizing
	std::vector<UINT> creaseCount; // The number of neighboring hard edges per vert
	UINT numVertices;
	UINT numUnpinned;

	Relaxer(
		short borderBehavior,
		short hardEdgeBehavior,
		short groupEdgeBehavior,
		std::vector<std::vector<size_t>> rawNeighbors,
		std::vector<std::vector<UCHAR>> rawHardEdges,
		const std::vector<UCHAR> &rawVertData
	);

	void edgeProject(
		const FLOAT basePoints[][NUM_COMPS],
		FLOAT smoothPoints[][NUM_COMPS]
	) const;

	void quickLaplacianSmooth(
		FLOAT verts2d[][NUM_COMPS],
		const FLOAT taubinBias=1.0
	) const;

};
