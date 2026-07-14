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
	//static MObject iterations;
	//static MObject smoothAlpha;
	static MObject wrinkleFreqVal;
	static MObject wrinkleAmpVal;
	//static MObject compressionThreshold;
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

	static std::vector<float> s_hostQInv;       // [numTris * 4]
	static std::vector<float> s_hostNormals;    // [numTris * 6]
	static std::vector<int>   s_hostVertIdx;    // [numTris * 3]
	static std::vector<float> s_hostRestCross;  // [numTris]
	static bool               s_topologyReady;

	// Per-frame BFS results — uploaded to GPU each frame
	static std::vector<float> s_phase;          // [numVerts]
	static std::vector<float> s_strainMask;     // [numVerts]
	static std::vector<float> s_vertexAmps;     // [numVerts]
	static std::vector<float> s_normals;        // [numVerts * 3]
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
	// Topology buffers — uploaded once from CPU meshTopology
	MAutoCLMem fTriQInvBuffer;
	MAutoCLMem fTriNormalsBuffer;
	MAutoCLMem fTriVertIdxBuffer;
	MAutoCLMem fTriRestCrossBuffer;
	bool fTopologyUploaded = false;

	// Per-vertex accumulator buffers — zeroed and scatter-written by strain kernel each frame
	MAutoCLMem fVertStrainAccumBuffer;
	MAutoCLMem fVertDirAccumBuffer;
	MAutoCLMem fVertAreaDeltaBuffer;

	// Per-frame CPU BFS results — uploaded each frame before displacement kernel
	MAutoCLMem fVertexNormalsBuffer;
	MAutoCLMem fWrinklePhaseBuffer;
	MAutoCLMem fStrainMaskBuffer;
	MAutoCLMem fVertexAmpsBuffer;

	// Helper methods
	bool extractLocatorMatrix(MDataBlock& block);

	// Pass 1: one thread per triangle — scatter strain into vertex accumulators
	cl_int enqueueStrainPass(
		MAutoCLEvent& syncEvent,
		const MGPUDeformerBuffer& inputPositions,
		unsigned int numTris,
		float warpStiffness, float weftStiffness, float areaStiffness,
		float frequency, float amplitude
	);

	// Pass 2: one thread per vertex — apply BFS phase + displacement
	cl_int enqueueDisplacePass(
		MAutoCLEvent& syncEvent,
		const MGPUDeformerBuffer& inputPositions,
		MGPUDeformerBuffer& outputPositions,
		float envelope, float threshold
	);

	// GPU data
	cl_float16        fLocatorMatrix;
	MOpenCLKernelInfo fStrainKernelInfo;   // strainKernel.cl
	MOpenCLKernelInfo fKernelInfo;         // kernel.cl (displacement)
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