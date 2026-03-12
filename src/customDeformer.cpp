#include "customDeformer.h"
#include "meshTopology.h"

#include <maya/MFnMatrixData.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>


MTypeId customDeformer::id(0x8000f);
MObject customDeformer::locatorMatrix;
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

	attributeAffects(customDeformer::locatorMatrix, customDeformer::outputGeom);
	
	return MStatus::kSuccess;
}

MStatus customDeformer::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& /*m*/, unsigned int multiIndex)
{
	MStatus returnStatus;

	float envelope = block.inputValue(customDeformer::envelope, &returnStatus).asFloat();
	if (envelope < 0.001f) return MStatus::kSuccess;

	int iterations = block.inputValue(customDeformer::iterations, &returnStatus).asInt();
	float smoothAlpha = block.inputValue(customDeformer::smoothAlpha, &returnStatus).asFloat();
	float wrinkleFreqVal = block.inputValue(customDeformer::wrinkleFreqVal, &returnStatus).asFloat();
	float wrinkleAmpVal = block.inputValue(customDeformer::wrinkleAmpVal, &returnStatus).asFloat();
	float compressionThreshold = block.inputValue(customDeformer::compressionThreshold, &returnStatus).asFloat();

	int numVerts = iter.count();

	MArrayDataHandle inputArray = block.inputArrayValue(MPxGeometryFilter::input, &returnStatus);
	returnStatus = inputArray.jumpToElement(multiIndex);
	MDataHandle inputHandle = inputArray.inputValue(&returnStatus);
	MObject currentInputMesh = inputHandle.child(MPxGeometryFilter::inputGeom).asMesh();

	MFnMesh currentMeshFn(currentInputMesh);


	if (!mInitialized) {
		mesh.buildFromMesh(currentInputMesh);
		mInitialized = true;
	}

	std::vector<MPoint> currentPos(numVerts);
	std::vector<MVector> currentNormals(numVerts);

	for (iter.reset(); !iter.isDone(); iter.next()) {
		int idx = iter.index();
		currentPos[idx] = iter.position();
		currentNormals[idx] = iter.normal();
	}

	std::vector<MPoint> readPos = currentPos;
	std::vector<MPoint> writePos = currentPos;

	/*
	  - Compute F = P * Q⁻¹  (deformation gradient); P  = [p1-p0 | p2-p0]
	  - Compute S = FᵀF       (strain tensor)
	  - Extract strain values  (S₁₁-1, S₂₂-1, S₁₂) as wrinkle magnitude/direction mask
	  - Displace vertices along normal, weighted by strain
	*/

	MItMeshPolygon polyIter(currentInputMesh, &returnStatus);

	MPointArray trianglePoints;
	MIntArray triangleVertexIndices;

	for (; !polyIter.isDone(); polyIter.next()) {
		polyIter.getTriangles(trianglePoints, triangleVertexIndices, MSpace::kObject);
	}

	for (int i = 0; i < mesh.tritoQInv.size(); ++i) {
		TriangleData triData = mesh.tritoQInv[i];

		// get the vertex ids for for the current tri
		int vertIdx0 = triData.vertIdx[0];
		int vertIdx1 = triData.vertIdx[1];
		int vertIdx2 = triData.vertIdx[2];

		// [ q00  q01
		//   q10  q11]
		float q00 = triData.qInv[0][0];
		float q01 = triData.qInv[0][1];
		float q10 = triData.qInv[1][0];
		float q11 = triData.qInv[1][1];

		MPoint aPos = trianglePoints[vertIdx0];
		MPoint bPos = trianglePoints[vertIdx1];
		MPoint cPos = trianglePoints[vertIdx2];

		//P = [p1 - p0 | p2 - p0]
		
		// [ p00  p01       [ q00  q01
		//   p10  p11   X     q10  q11 ]
		//   p20  p21]
		float p00 = bPos[0] - aPos[0];
		float p01 = cPos[0] - aPos[0];
		float p10 = bPos[1] - aPos[1];
		float p11 = cPos[1] - aPos[1];
		float p20 = bPos[2] - aPos[2];
		float p21 = cPos[2] - aPos[2];

		// Compute F = P * Q⁻¹(deformation gradient)
		float f00 = p00 * q00 + p01 * q10;
		float f01 = p00 * q01 + p01 * q11;
		float f10 = p10 * q00 + p11 * q10;
		float f11 = p10 * q01 + p11 * q11;
		float f20 = p20 * q00 + p21 * q10;
		float f21 = p20 * q01 + p21 * q11;


		// [ f00  f01			
		//   f10  f11   => F^T => [ f00  f10  f20
		//   f20  f21]				f01  f11  f21]	

		float t00 = f00;
		float t01 = f10;
		float t02 = f20;
		float t10 = f01;
		float t11 = f11;
		float t12 = f21;

		// S = FᵀF (strain tensor)	
		//							[ f00  f01
		// [ t00  t01  t02   X        f10  f11
		//   t10  t11  t12 ]	      f20  f21 ]

		float s00 = t00 * f00 + t01 * f10 + t02 * f20;
		float s01 = t00 * f01 + t01 * f11 + t02 * f21;
		float s10 = t10 * f00 + t11 * f11 + t12 * f20;
		float s11 = t10 * f01 + t11 * f11 + t12 * f21;

		//(S₁₁ - 1, S₂₂ - 1, S₁₂) as wrinkle magnitude / direction mask
		// S00 -1, S11 - 1, S01)
		MVector wrinkleMask(s00 - 1.0f, s11 - 1.0f, s01);

		MVector normal = currentNormals[i];

		float weight = normal * wrinkleMask;

		
		MVector tempPos = currentPos[i];
		writePos[i].x = tempPos.x + (normal * wrinkleMask * envelope);
		writePos[i].y = tempPos.y + (normal * wrinkleMask * envelope);
		writePos[i].z = tempPos.z + (normal * wrinkleMask * envelope);

	}

	iter.reset();
	for (; !iter.isDone(); iter.next()) {
		iter.setPosition(writePos[iter.index()]);
	}

	/*
	for (int iterCount = 0; iterCount < iterations; ++iterCount) {
	#pragma omp parallel for
		for (int i = 0; i < numVerts; ++i) {
			if (i >= mesh.vertIDConnections.size()) continue;

			MPoint myPos = readPos[i];
			MPoint neighborSum(0, 0, 0);
			float tensionSum = 0.0f;
			int valence = 0;

			for (auto neighborIdx : mesh.vertIDConnections[i]) {
				MPoint neighborPos = readPos[neighborIdx];
				neighborSum += neighborPos;

				double currentEdgeLen = myPos.distanceTo(neighborPos);
				std::pair<int, int> edgePair = std::make_pair(std::min(i, neighborIdx), std::max(i, neighborIdx));

				
				float restLen = mesh.vertstoRestLen[i];
				if (restLen > 0.0001f) {
					float normalizedTension = float(currentEdgeLen - restLen) / restLen;
					if (normalizedTension > 0) tensionSum += normalizedTension;
				}
	
				valence++;
			}

			if (valence > 0) {
				MPoint avgPos = neighborSum / valence;
				float avgTension = tensionSum / valence;
				MVector smoothVec = avgPos - myPos;

				float dynamicStiffness = avgTension * envelope * smoothAlpha;
				if (dynamicStiffness > 1.0f) dynamicStiffness = 1.0f;

				writePos[i] = myPos + (smoothVec * dynamicStiffness);
			}
		}
		readPos = writePos; // Swap buffer
	}


	#pragma omp parallel for
	for (int i = 0; i < numVerts; ++i) {
		if (i >= mesh.vertIDConnections.size()) continue;

		MPoint myPos = readPos[i];
		MVector myNormal = currentNormals[i];

		MVector W(1, 0, 0);

		MVector C = (W ^ myNormal).normal();

		// measure Anisotropic Compression
		float projectedCompression = 0.0f;
		float totalWeight = 0.0f;

		for (auto neighborIdx : mesh.vertIDConnections[i]) {
			MVector edgeDir = (readPos[neighborIdx] - myPos).normal();
			double currentEdgeLen = myPos.distanceTo(readPos[neighborIdx]);

			std::pair<int, int> edgePair = std::make_pair(std::min(i, neighborIdx), std::max(i, neighborIdx));
			
			float restLen = mesh.vertstoRestLen[i];
			if (restLen > 0.0001f) {
					
				float normalizedTension = float(currentEdgeLen - restLen) / restLen;
				float alignmentWeight = abs(edgeDir * C);

				if (normalizedTension < 0) {
					projectedCompression += (normalizedTension * alignmentWeight);
				}
				totalWeight += alignmentWeight;
			}
		
		}

		if (totalWeight > 0.0001f) {
			projectedCompression /= totalWeight;
		}

		MVector wrinkleOffset(1, 0, 0);

		if (projectedCompression < compressionThreshold) {
			float compressionFactor = (compressionThreshold - projectedCompression);

			MVector tempPos = myPos;
			double phase = (tempPos * C);
			double wave = sin(phase * wrinkleFreqVal);

			wrinkleOffset = myNormal * (wave * wrinkleAmpVal * compressionFactor * envelope);
		}

		writePos[i] = myPos + wrinkleOffset;
	}

	iter.reset();
	for (; !iter.isDone(); iter.next()) {
		iter.setPosition(writePos[iter.index()]);
	}
	*/
	return MStatus::kSuccess;
}