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
#include <maya/MPxDeformerNode.h> 
#include <vector>
#include "fastRelax.h"

#define DEFORMER_NAME "BlurRelax"
#define CHECKSTAT(m) if (!status) {status.perror(m); return status;};

// Double vs Float
typedef MPoint point_t;
typedef MPointArray pointArray_t;

void loadMayaTopologyData(
	MObject &mesh,
	MItGeometry& vertIter,
	std::vector<std::vector<size_t>> &neighbors, // A vector of neighbor indices per vertex
	std::vector<std::vector<char>> &hardEdges, // Bitwise per-neighbor data: edge is hard, edge along boundary
	std::vector<char> &vertData // Bitwise per-vert data: Group membership, geo boundary, group boundary,
);

class BlurRelax : public MPxDeformerNode {
	public:
		BlurRelax();
		~BlurRelax() override;

		static void* creator();
		static MStatus initialize();

		MStatus deform(MDataBlock& dataBlock, MItGeometry& vertIter, const MMatrix& matrix, UINT multiIndex) override;

	public:
		// local node attributes
		static MObject aIterations;
		static MObject aBorderBehavior;
		static MObject aHardEdgeBehavior;
		static MObject aGroupEdgeBehavior;
		static MObject aReproject;
		static MObject aTaubinBias;
		static MObject aRecomputeTopo;
		static MObject aDeltaMush;
		static MObject aDeltaMult;
		static MObject aDeltaBase;
		static MTypeId id;
	private:
		// Hash checking variables
		int hNumVerts = 0;
		int hNumPolys = 0;
		int hNumEdges = 0;
		short bbCheck = 255;
		short hbCheck = 255;
		short gbCheck = 255;

		// storage for this data doesn't change unless the hashes do
		std::vector<size_t> order;
		std::vector<size_t> invOrder;
		std::vector<std::vector<size_t>> neighbors; // A vector of neighbor indices per vertex
		std::vector<std::vector<char>> hardEdges; // Bitwise per-neighbor data: edge is hard, edge along boundary
		std::vector<char> vertData; // Bitwise per-vert data: Group membership, geo boundary, group boundary,
		std::vector<FLOAT> shiftVal; // normally 0.5; but it's 0.25 if on a hard edge
		std::vector<FLOAT> valence; // as float for vectorizing
		std::vector<UINT> creaseCount;

		MStatus getTrueWeights(
				MObject &mesh,
				MDataBlock &dataBlock,
				UINT index,
				std::vector<float> &weightVals,
				float envelope
				) const;

		void quickRelax(
			MObject &mesh,
			const short borderBehavior,
			const short hardEdgeBehavior,
			const short groupEdgeBehavior,
			const bool reproject,
			const float taubinBias,
			const FLOAT iterations,
			const UINT numVerts,
			const std::vector<char> &group,
			FLOAT(*verts)[4] // already resized
		);
};



