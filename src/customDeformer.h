#ifndef CUSTOM_DEFORMER
#define CUSTOM_DEFORMER

#include <maya/MPoint.h>
#include <maya/MGlobal.h>
#include <maya/MGeometry.h>

#include <maya/MPxDeformerNode.h>
#include <maya/MItGeometry.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include "meshTopology.h"

class customDeformer : public MPxDeformerNode {

public:
	customDeformer() {};
	virtual ~customDeformer() {};
	static void* creator();
	static MStatus initialize();
	virtual MStatus deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex);

	static MObject locatorMatrix;
	static MTypeId id;
	static MObject iterations;
	static MObject smoothAlpha;
	static MObject wrinkleFreqVal;
	static MObject wrinkleAmpVal;
	static MObject compressionThreshold;
	bool mInitialized = false;
	meshTopology mesh;

};


#endif
