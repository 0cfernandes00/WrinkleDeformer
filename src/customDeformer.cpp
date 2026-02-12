#include "customDeformer.h"
#include "meshTopology.h"

#include <maya/MFnMatrixData.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnNumericAttribute.h>


MTypeId customDeformer::id(0x8000f);
MObject customDeformer::locatorMatrix;
MObject customDeformer::angle;



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
	
	envelope = block.inputValue(customDeformer::envelope, &returnStatus).asFloat();
	

	MGlobal::displayInfo("printing inside deform block");

	MArrayDataHandle inputArray = block.inputArrayValue(MPxGeometryFilter::input, &returnStatus);
	returnStatus = inputArray.jumpToElement(multiIndex);
	MDataHandle inputData = inputArray.inputValue(&returnStatus);

	MDataHandle inputGeomHandle = inputData.child(MPxGeometryFilter::inputGeom);
	MObject inputMesh = inputGeomHandle.asMesh();

	std::vector<MPoint> currentPos;
	currentPos.clear();

	std::map<std::pair<int, int>, float> tension;


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

		MGlobal::displayInfo("INITIALIZING FROM ORIGINAL MESH");
		mInitialized = true;

	}

	int numVerts = iter.count();
	if (currentPos.size() != numVerts) {
		currentPos.resize(numVerts);
	}

	// Copy all current positions to our vector so we can random-access neighbors later
	for (iter.reset(); !iter.isDone(); iter.next()) {
		currentPos[iter.index()] = iter.position();
	}

	iter.reset();

	for(; !iter.isDone(); iter.next()) {

		// iterate over connected verts map
		int currVertIdx = iter.index();

		if (currVertIdx >= mesh.vertIDConnections.size()) continue;

		std::vector<int>& neighbors = mesh.vertIDConnections[currVertIdx];
		float tensionSum = 0.0f;
		int valence = 0;

		MPoint myCurrentPos = currentPos[currVertIdx];

		for (int i = 0; i < neighbors.size(); ++i) {

			int vert = neighbors[i];

			// check it's not getting the same vert
			if (vert != currVertIdx) {
				
				// retrieve the tension
				float vertDist = myCurrentPos.distanceTo(currentPos[vert]);

				std::pair<int,int> edgePair = std::make_pair(std::min(currVertIdx, vert), std::max(currVertIdx, vert));

				if (mesh.vertstoRestLen.find(edgePair) != mesh.vertstoRestLen.end()) {
					float vertRL = mesh.vertstoRestLen[edgePair];

					if (vertRL > 0.0001) {
						tensionSum += float((vertDist - vertRL) / vertRL);
						valence++;
					}

				}
			}
		} 

		float avgTension = (valence > 0) ? (tensionSum / valence) : 0.0f;

		
		if (abs(avgTension) > 0.001f) {
			float outScalar = avgTension * envelope;

			MVector normal = iter.normal();

			MVector newPos = myCurrentPos + (normal * outScalar);

			// Positive = stretch, Negative = compression
			// Move Vertex
			// Stretch (pos tension) -> Pushes OUT (along normal)
			// Compress (neg tension) -> Pushes IN (inverse normal)
			iter.setPosition(newPos);
			
		}

	}


	return returnStatus;
}