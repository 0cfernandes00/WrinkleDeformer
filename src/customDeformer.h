#ifndef CUSTOM_DEFORMER
#define CUSTOM_DEFORMER

#include <maya/MPoint.h>
#include <maya/MGlobal.h>
#include <maya/MGeometry.h>
#include <maya/MPxGPUStandardDeformer.h>
#include <maya/MOpenCLUtils.h>
#include <maya/MGPUDeformerRegistry.h>

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
	static MObject warpStiffness;
	static MObject weftStiffness;
	static MObject areaStiffness;
	static MObject restMesh;
	bool mInitialized = false;
	meshTopology mesh;

	std::vector<std::vector<MVector>> m_threadAccum;
	std::vector<std::vector<int>> m_threadCount;
	std::vector<std::vector<MVector>> m_threadWrinkleDir;
	std::vector<std::vector<float>> m_threadPhysAmp;
	std::vector<std::vector<MVector>> m_threadAreaAccum; 

	std::vector<double> m_wrinklePhase;
	std::vector<float> m_strainMask;
	std::vector<MVector> m_vertexDirs;
	std::vector<float> m_vertexAmps;
	std::vector<MPoint> m_writePos;
	std::vector<float> m_pts;
	std::vector<float> m_nrms;
	std::vector<MPoint> m_currentPos;
	std::vector<MVector> m_vertexAreaDelta;

	std::vector<bool> m_visited;
	std::vector<bool> m_isBoundary;
	std::vector<int> m_currentFrontier;
	std::vector<int> m_nextFrontier;


	int m_cachedNumVerts = -1;
	int m_cachedNumThreads = -1;

	static MString pluginPath;

};

class wrinkleGPUDeformer : public MPxGPUStandardDeformer {

public:
	wrinkleGPUDeformer();
	virtual ~wrinkleGPUDeformer();

	MPxGPUDeformer::DeformerStatus evaluate(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& outputPlug, const MPlugArray& inputPlugs, const MGPUDeformerData& inputData, MGPUDeformerData& outputData) override;
	void terminate() override;
	static MGPUDeformerRegistrationInfo* getGPUDeformerInfo();
	static bool validateNodeInGraph(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, MStringArray* messages);
	static bool validateNodeValues(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, MStringArray* messages);

private:
	MAutoCLMem fTriQInvBuffer;
	MAutoCLMem fTriNormalsBuffer;
	MAutoCLMem fTriVertIdxBuffer;
	MAutoCLMem fTriRestCrossBuffer;
	MAutoCLMem fVertexNormalsBuffer;
	MAutoCLMem fWrinklePhaseBuffer;   
	MAutoCLMem fStrainMaskBuffer;    
	MAutoCLMem fVertexAmpsBuffer;    
	bool fTopologyUploaded = false;

	MOpenCLKernelInfo fStrainKernelInfo;
	MAutoCLMem fVertStrainAccumBuffer;  // GPU-side per-vertex accumulators
	MAutoCLMem fVertDirAccumBuffer;
	MAutoCLMem fVertAreaDeltaBuffer;


	// helper methods
	bool extractLocatorMatrix(MDataBlock& block);
	cl_int enqueueDeformation(
		MAutoCLEvent& syncEvent,
		const MGPUDeformerBuffer& inputPositions, const MGPUDeformerBuffer& triQInv, const MGPUDeformerBuffer& triNormals,
		const MGPUDeformerBuffer& triVertIdx, const MGPUDeformerBuffer& triRestCrossSqLen,
		MGPUDeformerBuffer& outputPositions, MGPUDeformerBuffer& vertStrainAccum, MGPUDeformerBuffer& vertDirAccum, MGPUDeformerBuffer& vertAreaDelta,
		float threshold, float amplitude, float frequency,
		float warpStiffness, float weftStiffness, float areaStiffness
	);	
	// Storage for data on the GPU
	cl_float16 fLocatorMatrix;
	// Kernel
	MOpenCLKernelInfo fKernelInfo;

};


class wrinkleDeformerNodeGPUDeformerInfo : public MGPUDeformerRegistrationInfo
{
public:
	wrinkleDeformerNodeGPUDeformerInfo();
	~wrinkleDeformerNodeGPUDeformerInfo() override {}
	MPxGPUStandardDeformer* createGPUDeformer() override
	{
		return new wrinkleGPUDeformer();
	}
	
	bool validateNodeInGraph(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages) override
	{
		return wrinkleGPUDeformer::validateNodeInGraph(block, evaluationNode, plug, messages);
	}
	bool validateNodeValues(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages) override
	{
		return wrinkleGPUDeformer::validateNodeValues(block, evaluationNode, plug, messages);
	}
	
};


#endif
