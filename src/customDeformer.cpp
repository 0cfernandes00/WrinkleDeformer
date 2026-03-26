#include "customDeformer.h"
#include "meshTopology.h"
#include <omp.h>
#include <queue>

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

	wrinkleFreqVal = nAttr.create("wrinkleFreq", "wrinkleFreq", MFnNumericData::kFloat, 5.0f);
	nAttr.setStorable(true);
	nAttr.setConnectable(true);
	nAttr.setKeyable(true);
	addAttribute(wrinkleFreqVal);

	wrinkleAmpVal = nAttr.create("wrinkleAmp", "wrinkleAmp", MFnNumericData::kFloat, 1.5f);
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
		mesh.buildFromMesh(currentInputMesh, numVerts);
		mInitialized = true;
	}

	//std::vector<MPoint> currentPos(numVerts);
	//std::vector<MVector> currentNormals(numVerts);
	/*
	for (iter.reset(); !iter.isDone(); iter.next()) {
		int idx = iter.index();
		//currentPos[idx] = iter.position();
		currentNormals[idx] = iter.normal();
	}
	*/


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

	if (numVerts != m_cachedNumVerts || numThreads != m_cachedNumThreads) {
		m_threadAccum.assign(numThreads, std::vector<MVector>(numVerts));
		m_threadCount.assign(numThreads, std::vector<int>(numVerts, 0));
		m_threadWrinkleDir.assign(numThreads, std::vector<MVector>(numVerts));
		m_threadPhysAmp.assign(numThreads, std::vector<float>(numVerts, 0.0f));

		m_wrinklePhase.assign(numVerts, -1.0);
		m_strainMask.assign(numVerts, 0.0f);
		m_vertexDirs.assign(numVerts, MVector(0, 0, 0));
		m_vertexAmps.assign(numVerts, 0.0f);

		m_cachedNumVerts = numVerts;
		m_cachedNumThreads = numThreads;

		m_pts.resize(numVerts * 3);
		m_nrms.resize(numVerts * 3);
		m_currentPos.resize(numVerts);
		m_writePos.resize(numVerts);

		m_visited.assign(numVerts, false);
		m_isBoundary.assign(numVerts, false);
		m_currentFrontier.reserve(numVerts);
		m_nextFrontier.reserve(numVerts);
	}

	for (int t = 0; t < numThreads; ++t) {
		std::fill(m_threadAccum[t].begin(), m_threadAccum[t].end(), MVector(0, 0, 0));
		std::fill(m_threadCount[t].begin(), m_threadCount[t].end(), 0);
		std::fill(m_threadWrinkleDir[t].begin(), m_threadWrinkleDir[t].end(), MVector(0, 0, 0));
		std::fill(m_threadPhysAmp[t].begin(), m_threadPhysAmp[t].end(), 0.0f);
	}

	std::fill(m_wrinklePhase.begin(), m_wrinklePhase.end(), -1.0);
	std::fill(m_strainMask.begin(), m_strainMask.end(), 0.0f);
	std::fill(m_vertexDirs.begin(), m_vertexDirs.end(), MVector(0, 0, 0));
	std::fill(m_vertexAmps.begin(), m_vertexAmps.end(), 0.0f);
	std::fill(m_visited.begin(), m_visited.end(), false);
	std::fill(m_isBoundary.begin(), m_isBoundary.end(), false);
	m_currentFrontier.clear();
	m_nextFrontier.clear();

	m_writePos = m_currentPos;

	for (int k = 0; k < numVerts; ++k) {
		m_pts[k * 3] = allPoints[k].x;
		m_pts[k * 3 + 1] = allPoints[k].y;
		m_pts[k * 3 + 2] = allPoints[k].z;
		m_nrms[k * 3] = allNormals[k].x;
		m_nrms[k * 3 + 1] = allNormals[k].y;
		m_nrms[k * 3 + 2] = allNormals[k].z;
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numTris; ++i) {
		int tid = omp_get_thread_num();
		const TriangleData& triData = mesh.tritoQInv[i];

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

		MPoint aPos = MPoint(m_pts[vertIdx0*3],m_pts[vertIdx0*3 +1], m_pts[vertIdx0*3+2]);
		MPoint bPos = MPoint(m_pts[vertIdx1*3],m_pts[vertIdx1*3+1], m_pts[vertIdx1*3+2]);
		MPoint cPos = MPoint(m_pts[vertIdx2*3],m_pts[vertIdx2*3+1], m_pts[vertIdx2*3+2]);

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
		// Principal direction via simple 2x2 eigenvector
		float angle = 0.5f * atan2(2.0f * s01, s00 - s11);

		// Get minimum principal strain - the most compressed direction
		// From the 2x2 eigenvalue decomposition of S
		float trace = s00 + s11;


		float det = s00 * s11 - s01 * s10;
		if (fabs(det) < 1e-6f) continue;

		float disc = sqrt(std::max(0.0f, (trace * trace * 0.25f) - det));
		float eigenMin = (trace * 0.5f) - disc; // most compressed eigenvalue

		// Compression factor - how much the surface shrank in that direction
		float compressionFactor = sqrt(std::max(0.0f, eigenMin));
		float sMag = 1.0f - compressionFactor;

		// Only apply if it's actually compressed
		strainMagnitude = (sMag > 0.0f) ? sMag : 0.0f;

		float waveLen = 1.0f / wrinkleFreqVal;
		float physicalAmplitude = 0.0f;
		if (compressionFactor < 0.0001f || compressionFactor >= 1.0f) {
			physicalAmplitude = 0.0f;
		}
		else {
			float invC = 1.0f / compressionFactor;
			physicalAmplitude = (waveLen / (2.0f * M_PI)) *
				sqrt(std::max(0.0f, invC * invC - 1.0f));
		}

		physicalAmplitude *= wrinkleAmpVal;
	
		float compX = cos(angle);
		float compY = sin(angle);

		float wrinkleX = -compY;
		float wrinkleY = compX;

		// Transform wrinkle direction from UV space to world space using n1, n2
		MVector n1 = triData.normal[0];
		MVector n2 = triData.normal[1];
		MVector wrinkleDirWorld = (n1 * wrinkleX + n2 * wrinkleY);
		double dirLen = wrinkleDirWorld.length();
		if (dirLen < 1e-6) continue;
		wrinkleDirWorld = wrinkleDirWorld / dirLen;
		//wrinkleDirWorld = MVector(1, 0, 0);

		//strainMagnitude *= triData.windingSign;

		//if (strainMagnitude >= 0.0f) continue;

		for (int j = 0; j < 3; ++j) {
		
			int vertIdx = triData.vertIdx[j];
			if (vertIdx >= numVerts) continue;

			if (m_threadCount[tid][vertIdx] > 0) {
				// Flip if pointing roughly opposite to what's already accumulated
				MVector existing = m_threadWrinkleDir[tid][vertIdx];
				if ((existing * wrinkleDirWorld) < 0.0f) {
					wrinkleDirWorld = -wrinkleDirWorld;
				}
			}

			m_threadAccum[tid][vertIdx] += MVector(m_nrms[vertIdx*3], m_nrms[vertIdx*3+1], m_nrms[vertIdx*3+2]) * strainMagnitude;
			m_threadWrinkleDir[tid][vertIdx] += wrinkleDirWorld;
			m_threadCount[tid][vertIdx]++;
			m_threadPhysAmp[tid][vertIdx] += physicalAmplitude;
		}

	}
	
	#pragma omp parallel for schedule(static)
	for (int k = 0; k < numVerts; ++k) {
		MVector accum(0, 0, 0);
		MVector dirAccum(0, 0, 0);
		float physicalAmplitude = 0.0f;
		int count = 0;
		for (int t = 0; t < numThreads; ++t) {
			accum += m_threadAccum[t][k];
			count += m_threadCount[t][k];
			dirAccum += m_threadWrinkleDir[t][k];
			physicalAmplitude += m_threadPhysAmp[t][k];
		}
		if (count > 0) {
			accum /= count;
			dirAccum = dirAccum.normal();
			m_strainMask[k] = accum.length();
			m_vertexDirs[k] = dirAccum;
			m_vertexAmps[k] = physicalAmplitude / count;
		}
	}

	/*	BFS Search	*/
	// 1. Clear frontiers and reset visited state
	m_currentFrontier.clear();
	std::fill(m_visited.begin(), m_visited.end(), false);
	std::fill(m_isBoundary.begin(), m_isBoundary.end(), false);

	// 2. Unified Seeding Logic
	for (int k = 0; k < numVerts; ++k) {
		// Look for "uncompressed" vertices (using your threshold, not a hardcoded 0.0001f)
		if (m_strainMask[k] <= 0.0001f) {
			int start = mesh.adjacencyStart[k];
			int count = mesh.adjacencyCount[k];

			for (int n = 0; n < count; ++n) {
				int neighborIdx = mesh.adjacencyData[start + n];

				// If an uncompressed vertex touches a compressed one, it's a seed!
				if (m_strainMask[neighborIdx] > 0.0001f) {
					m_isBoundary[k] = true;
					m_wrinklePhase[k] = 0.0; // Initialize phase
					m_currentFrontier.push_back(k);
					m_visited[k] = true;
					break;
				}
			}
		}
	}

	std::vector<int> hopsFromBoundary(numVerts, -1);
	// seed
	for (int k = 0; k < numVerts; ++k)
		if (m_isBoundary[k]) hopsFromBoundary[k] = 0;

	// BFS counting hops
	// then print
	int maxHops = 0;
	for (int k = 0; k < numVerts; ++k)
		if (m_strainMask[k] > 0.0001f && hopsFromBoundary[k] > maxHops)
			maxHops = hopsFromBoundary[k];
	MGlobal::displayInfo(MString("Max hops into compressed region: ") + maxHops);

	int compressedCount = 0;
	int boundaryCount = 0;
	float maxAmp = 0.0f;

	for (int k = 0; k < numVerts; ++k) {
		if (m_strainMask[k] > 0.0001f) compressedCount++;
		if (m_vertexAmps[k] > maxAmp) maxAmp = m_vertexAmps[k];
	}
	// after BFS
	for (int k = 0; k < numVerts; ++k) {
		if (m_isBoundary[k]) boundaryCount++;
	}

	MGlobal::displayInfo(
		MString("Compressed: ") + compressedCount +
		MString(" Boundary: ") + boundaryCount +
		MString(" MaxAmp: ") + maxAmp
	);

	// Propagate phase through connectivity
	while (!m_currentFrontier.empty()) {
		m_nextFrontier.clear();
	
		for (int current : m_currentFrontier) {
			MPoint currentPos = MPoint(m_pts[current * 3], m_pts[current * 3 + 1], m_pts[current * 3 + 2]);
			MVector currentDir = m_vertexDirs[current];

			int start = mesh.adjacencyStart[current];
			int count = mesh.adjacencyCount[current];
			for (int n = 0; n < count; ++n) {
				int neighborIdx = mesh.adjacencyData[start + n];
				
				if (m_visited[neighborIdx]) continue;
				m_visited[neighborIdx] = true;

				MPoint neighborPos = MPoint(m_pts[neighborIdx * 3], m_pts[neighborIdx * 3 + 1], m_pts[neighborIdx * 3 + 2]);
				MVector edge = MVector(neighborPos - currentPos);

				//MVector tmp = m_vertexDirs[neighborIdx];
				//tmp = normal * strainMask[k];
				MVector dirToUse = currentDir;
				if (dirToUse.length() < 0.001) {
					dirToUse = MVector(1, 0, 0); // Fallback
				}
				//MVector dirToUse = (currentDir.length() > 0.0001f) ?
					//currentDir : m_vertexDirs[neighborIdx];
				dirToUse.normalize();

				// Project edge onto wrinkle direction to get phase increment
				double phaseIncrement = (edge * dirToUse) * wrinkleFreqVal;
				m_wrinklePhase[neighborIdx] = m_wrinklePhase[current] + phaseIncrement;
				m_nextFrontier.push_back(neighborIdx);
			}
		}
		std::swap(m_currentFrontier, m_nextFrontier);
	}

	for (int k = 0; k < numVerts; ++k) {

		if (m_strainMask[k] > 0.0001f && m_visited[k]) {
			MVector normal = MVector(m_nrms[k*3], m_nrms[k*3+1], m_nrms[k*3+2]);
			MPoint pos = MPoint(m_pts[k*3], m_pts[k*3 + 1], m_pts[k*3 + 2]);
			//MVector tempPos = MVector(pos.x, pos.y, pos.z);

			double phase = m_wrinklePhase[k];
			double wave = pow(abs(sin(phase)), 0.5) * (sin(phase) > 0 ? 1.0 : -1.0);

			float wrinkleDisp = (float)(wave * m_vertexAmps[k] * envelope);
			m_writePos[k] = pos + normal * wrinkleDisp;
		}
		else {
			m_writePos[k] = MPoint(m_pts[k * 3], m_pts[k * 3 + 1], m_pts[k * 3 + 2]);
		}
	}
	
	iter.reset();
	for (; !iter.isDone(); iter.next()) {
		iter.setPosition(m_writePos[iter.index()]);
	}
	return MStatus::kSuccess;
}