#include <meshTopology.h>


void meshTopology::buildFromMesh(MObject& mesh) {


	MStatus status;
	MItMeshEdge edgeIter(mesh, &status);


	//std::map<std::pair<MPoint*, MPoint*>, float> vertstoRestLen;
	for (; !edgeIter.isDone(); edgeIter.next()) {
		// get vertex indices for the edge
		int vertId1 = edgeIter.index(0);
		int vertId2 = edgeIter.index(1);

		// store connectivity for each vertex
		vertIDConnections[vertId1].push_back(vertId2);
		vertIDConnections[vertId2].push_back(vertId1);

		MPoint pointId1 = edgeIter.point(0, MSpace::kObject, &status);
		MPoint pointId2 = edgeIter.point(1, MSpace::kObject, &status);

		double dist = pointId1.distanceTo(pointId2);

		std::pair pointPair = std::make_pair(std::min(vertId1, vertId2), std::max(vertId1, vertId2));
		vertstoRestLen[pointPair] = dist;
	}

	MGlobal::displayInfo(MString("Stored ") + vertstoRestLen.size() + MString(" edges"));
	if (vertstoRestLen.size() > 0) {
		auto firstEdge = vertstoRestLen.begin();
		MGlobal::displayInfo(MString("First edge rest length: ") + firstEdge->second);
	}

}

/*
const Vertex& meshTopology::getVertex(uint32_t index) const {

	Vertex vert;
	return vert;
}
*/