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

#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include "blurRelaxNode.h"
#include "version.h"

MStatus initializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj, "Blur Studio", VERSION_STRING, "Any");
	result = plugin.registerNode(DEFORMER_NAME, BlurRelax::id, BlurRelax::creator, BlurRelax::initialize, MPxNode::kDeformerNode);
	if (MGlobal::mayaState() == MGlobal::kInteractive) {
		MGlobal::executeCommand("makePaintable -attrType \"multiFloat\" -sm \"deformer\" \"" DEFORMER_NAME "\" \"weights\";");
		MGlobal::executePythonCommandOnIdle("import BlurRelax.menu");
		MGlobal::executePythonCommandOnIdle("BlurRelax.menu.create_menuitems()");
	}
	return result;
}

MStatus uninitializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj);
	result = plugin.deregisterNode(BlurRelax::id);
	return result;
}
