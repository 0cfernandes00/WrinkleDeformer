#include "customDeformer.h"
#include "meshTopology.h"
#include <omp.h>

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

	//std::vector<MPoint> readPos = currentPos;
	std::vector<MPoint> writePos = currentPos;

	/*
	  - Compute F = P * Q⁻¹  (deformation gradient); P  = [p1-p0 | p2-p0]
	  - Compute S = FᵀF       (strain tensor)
	  - Extract strain values  (S₁₁-1, S₂₂-1, S₁₂) as wrinkle magnitude/direction mask
	  - Displace vertices along normal, weighted by strain
	*/

	MPointArray allPoints;
	currentMeshFn.getPoints(allPoints, MSpace::kObject);

	MFloatVectorArray allNormals;
	currentMeshFn.getVertexNormals(false, allNormals, MSpace::kObject);

	int numTris = mesh.tritoQInv.size();

	int numThreads = omp_get_max_threads();
	std::vector<std::vector<MVector>> threadAccum(numThreads,
		std::vector<MVector>(numVerts, MVector(0, 0, 0)));
	std::vector<std::vector<int>> threadCount(numThreads,
		std::vector<int>(numVerts, 0));


	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numTris; ++i) {
		int tid = omp_get_thread_num();
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

		MPoint aPos = allPoints[vertIdx0];
		MPoint bPos = allPoints[vertIdx1];
		MPoint cPos = allPoints[vertIdx2];

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
		float s10 = t10 * f00 + t11 * f10 + t12 * f20;
		float s11 = t10 * f01 + t11 * f11 + t12 * f21;

		//(S₁₁ - 1, S₂₂ - 1, S₁₂) as wrinkle magnitude / direction mask
		// (S00 -1, S11 - 1, S01)
		float strainMagnitude = (s00 - 1.0f) + (s11 - 1.0f);

		if (strainMagnitude >= 0.0f) continue;

		for (int j = 0; j < 3; ++j) {
		
			int vertIdx = triData.vertIdx[j];
			if (vertIdx > numVerts) continue;

			threadAccum[tid][vertIdx] += MVector(allNormals[vertIdx]) * strainMagnitude;
			threadCount[tid][vertIdx]++;
		}
	}
	
	for (int k = 0; k < numVerts; ++k) {
		MVector accum(0, 0, 0);
		int count = 0;
		for (int t = 0; t < numThreads; ++t) {
			accum += threadAccum[t][k];
			count += threadCount[t][k];
		}
		if (count > 0) {
			accum /= count;
		}
		float strainWeight = accum.length(); // scalar compression magnitude

		if (strainWeight > 0.0001f) {
			MVector normal = MVector(allNormals[k]);

			// Use strain magnitude to gate the wrinkle
			// Use position projected onto a wrinkle direction for the wave phase
			MPoint pos = allPoints[k];

			// Wrinkle direction - across the compression axis
			// For now use world X as wrinkle direction, can be made anisotropic later
			double phase = pos.x * wrinkleFreqVal;
			double wave = sin(phase);

			float wrinkleDisp = (float)(wave * wrinkleAmpVal * strainWeight * envelope);
			writePos[k] = allPoints[k] + normal * wrinkleDisp;
		}
		else {
			writePos[k] = allPoints[k];
		}
	}

	int compressedVerts = 0;
	float maxDisp = 0.0f;
	for (int k = 0; k < numVerts; ++k) {
		MVector d = writePos[k] - allPoints[k];
		float len = d.length();
		if (len > 0.0001f) compressedVerts++;
		if (len > maxDisp) maxDisp = len;
	}
	MGlobal::displayInfo(MString("Compressed verts: ") + compressedVerts +
		MString(" MaxDisp: ") + maxDisp);

	iter.reset();
	for (; !iter.isDone(); iter.next()) {
		iter.setPosition(writePos[iter.index()]);
	}
	return MStatus::kSuccess;
}