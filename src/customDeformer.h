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

	std::vector<std::vector<MVector>> m_threadAccum;
	std::vector<std::vector<int>> m_threadCount;
	std::vector<std::vector<MVector>> m_threadWrinkleDir;
	std::vector<std::vector<float>> m_threadPhysAmp;

	std::vector<double> m_wrinklePhase;
	std::vector<float> m_strainMask;
	std::vector<MVector> m_vertexDirs;
	std::vector<float> m_vertexAmps;
	std::vector<MPoint> m_writePos;
	std::vector<float> m_pts;
	std::vector<float> m_nrms;
	std::vector<MPoint> m_currentPos;

	std::vector<bool> m_visited;
	std::vector<bool> m_isBoundary;
	std::vector<int> m_currentFrontier;
	std::vector<int> m_nextFrontier;


	int m_cachedNumVerts = -1;
	int m_cachedNumThreads = -1;

};


#endif
