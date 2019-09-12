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
#include "fastMayaRelax.h"

#define DEFORMER_NAME "BlurRelax"
#define CHECKSTAT(m) if (!status) {status.perror(m); return status;};

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

		MayaRelaxer *relaxer;

		MStatus getTrueWeights(
			MObject &mesh, MDataBlock &dataBlock, UINT index,
			std::vector<float> &weightVals, float envelope) const;
};

