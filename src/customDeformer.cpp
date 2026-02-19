#include "customDeformer.h"
#include "meshTopology.h"

#include <maya/MFnMatrixData.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>


MTypeId customDeformer::id(0x8000f);
MObject customDeformer::locatorMatrix;
MObject customDeformer::angle;
MObject customDeformer::iterations;
MObject customDeformer::smoothAlpha;
MObject customDeformer::wrinkleFreqVal;
MObject customDeformer::wrinkleAmpVal;
MObject customDeformer::compressionThreshold;



void* customDeformer::creator() 
{
	return new customDeformer();
}

MStatus customDeformer::initialize() 
{
	// define attrs
	MFnMatrixAttribute mAttr;
	MFnNumericAttribute nAttr;

	locatorMatrix = mAttr.create("locatorMatrix", "lm");
	mAttr.setStorable(false);
	mAttr.setConnectable(true);
	addAttribute(locatorMatrix);

	angle = nAttr.create("angle", "ang", MFnNumericData::kFloat, 0);
	nAttr.setStorable(true);
	nAttr.setConnectable(true);
	nAttr.setKeyable(true);
	addAttribute(angle);

	iterations = nAttr.create("iterations", "iterations", MFnNumericData::kInt, 10);
	nAttr.setStorable(true);
	nAttr.setConnectable(true);
	nAttr.setKeyable(true);
	addAttribute(iterations);

	smoothAlpha = nAttr.create("smoothAlpha", "sa", MFnNumericData::kFloat, 0.1f);
	nAttr.setStorable(true);
	nAttr.setConnectable(true);
	nAttr.setKeyable(true);
	addAttribute(smoothAlpha);

	wrinkleFreqVal = nAttr.create("wrinkleFreq", "wrinkleFreq", MFnNumericData::kFloat, 10.0f);
	nAttr.setStorable(true);
	nAttr.setConnectable(true);
	nAttr.setKeyable(true);
	addAttribute(wrinkleFreqVal);

	wrinkleAmpVal = nAttr.create("wrinkleAmp", "wrinkleAmp", MFnNumericData::kFloat, 0.5f);
	nAttr.setStorable(true);
	nAttr.setConnectable(true);
	nAttr.setKeyable(true);
	addAttribute(wrinkleAmpVal);

	compressionThreshold = nAttr.create("compressionThreshold", "compressionThreshold", MFnNumericData::kFloat, -0.1f);
	nAttr.setStorable(true);
	nAttr.setConnectable(true);
	nAttr.setKeyable(true);
	addAttribute(compressionThreshold);

	
	attributeAffects(customDeformer::iterations, customDeformer::outputGeom);
	attributeAffects(customDeformer::smoothAlpha, customDeformer::outputGeom);
	attributeAffects(customDeformer::wrinkleFreqVal, customDeformer::outputGeom);
	attributeAffects(customDeformer::wrinkleAmpVal, customDeformer::outputGeom);
	attributeAffects(customDeformer::compressionThreshold, customDeformer::outputGeom);

	attributeAffects(customDeformer::angle, customDeformer::outputGeom);
	attributeAffects(customDeformer::locatorMatrix, customDeformer::outputGeom);
	
	return MStatus::kSuccess;
}

MStatus customDeformer::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& /*m*/, unsigned int multiIndex)
{
	MStatus returnStatus;
	MDataHandle envData;
	float angle;
	float envelope;
	float iterations;
	float smoothAlpha;
	float wrinkleFreqVal;
	float wrinkleAmpVal;
	float compressionThreshold;
	
	envelope = block.inputValue(customDeformer::envelope, &returnStatus).asFloat();
	iterations = block.inputValue(customDeformer::iterations, &returnStatus).asInt();
	smoothAlpha = block.inputValue(customDeformer::smoothAlpha, &returnStatus).asFloat();
	wrinkleFreqVal = block.inputValue(customDeformer::wrinkleFreqVal, &returnStatus).asFloat();
	wrinkleAmpVal = block.inputValue(customDeformer::wrinkleAmpVal, &returnStatus).asFloat();
	compressionThreshold = block.inputValue(customDeformer::compressionThreshold, &returnStatus).asFloat();

	std::vector<MPoint> currentPos;
	std::vector<MVector> currentNormals;
	currentPos.clear();

	std::map<std::pair<int, int>, float> tension;

	if (envelope < 0.001f) return MStatus::kSuccess;

	int numVerts = iter.count();


	// store bind pose data
	if (!mInitialized) {

		// construct connectivity data (vertex pairs to each edge)
		MArrayDataHandle inputArray = block.inputArrayValue(MPxGeometryFilter::input, &returnStatus);
		returnStatus = inputArray.jumpToElement(multiIndex);
		MDataHandle inputHandle = inputArray.inputValue(&returnStatus);

		MDataHandle geomHandle = inputHandle.child(MPxGeometryFilter::inputGeom);
		MObject origMesh = geomHandle.asMesh();

		MFnMesh meshFn(origMesh);

		// Create a copy of the mesh data
		MObject restMeshData = meshFn.copy(origMesh);

		mesh.buildFromMesh(restMeshData);

		mInitialized = true;

	}

	if (currentPos.size() != numVerts || currentNormals.size() != numVerts) {
		currentPos.resize(numVerts);
		currentNormals.resize(numVerts);
	}

	for (iter.reset(); !iter.isDone(); iter.next()) {
		currentPos[iter.index()] = iter.position();
		currentNormals[iter.index()] = iter.normal();
	}

	
	// ping pong buffers smoothing
	std::vector<MPoint> readPos = currentPos;
	std::vector<MPoint> writePos = currentPos;

	for (int iterCount = 0; iterCount < iterations; ++iterCount) {

		for (int i = 0; i < numVerts; ++i) {

			if (i >= mesh.vertIDConnections.size()) continue;

			std::vector<int>& neighbors = mesh.vertIDConnections[i];

			MPoint myPos = readPos[i];
			MVector myNormal = currentNormals[i];

			MPoint neighborSum(0, 0, 0);
			float tensionSum = 0.0f;
			int valence = 0;


			MVector maxCompressionDir(0, 1, 0);
			float minTension = 0.0f;

			for (auto neighborIdx : neighbors) {

				MPoint neighborPos = readPos[neighborIdx];
				neighborSum += neighborPos;

				double currentEdgeLen = myPos.distanceTo(neighborPos);

				std::pair<int, int> edgePair = std::make_pair(std::min(i, neighborIdx), std::max(i, neighborIdx));
				
				if (mesh.vertstoRestLen.count(edgePair)) {
					float restLen = mesh.vertstoRestLen[edgePair];
					if (restLen > 0.0001f) {
						float edgeTension = float(currentEdgeLen - restLen);
						float normalizedTension = edgeTension / restLen;

						if (edgeTension > 0) {
							tensionSum += edgeTension;
						}
						if (normalizedTension < minTension) {
							minTension = normalizedTension;
							maxCompressionDir = (neighborPos - myPos).normal();
						}

					}
				}
				
				valence++;
			}
			
			MVector totalOffset(0, 0, 0);
			if (valence > 0) {
				MPoint avgPos = neighborSum / valence;
				float avgTension = tensionSum / valence;

				MVector smoothVec = avgPos - myPos;

				// only smooth stretched area

				float dynamicStiffness = avgTension * envelope;

				if (dynamicStiffness > 1.0f) dynamicStiffness = 1.0f;

				totalOffset += (smoothVec * dynamicStiffness);
			}

			float threshold = -0.1f; // User attribute: compressionThreshold
			if (minTension < threshold) {

				// A. Magnitude: How deep is the wrinkle?
				// Map tension (e.g. -0.5) to a positive multiplier (0.0 to 1.0)
				float compressionFactor = (compressionThreshold - minTension); // The more negative, the bigger this gets

				// B. Phase: Where are we on the wave?
				// We project our position onto the compression direction.
				// This ensures vertices "in a line" along the compression get different sine values.
				MVector finalPos = myPos;
				double phase = (finalPos * maxCompressionDir); // Dot Product acts as 1D projection

				// C. The Wave Function
				// sin(Position * Frequency) * Amplitude * CompressionStrength
				double wave = sin(phase * wrinkleFreqVal);

				// D. Displacement
				MVector displacement = myNormal * (wave * wrinkleAmpVal * compressionFactor * envelope);

				totalOffset += displacement;
			}

			writePos[i] = myPos + totalOffset;
		}
		readPos = writePos;
	}
	iter.reset();
	for (; !iter.isDone(); iter.next()) {
		iter.setPosition(writePos[iter.index()]);
	}
	

	return returnStatus;
}