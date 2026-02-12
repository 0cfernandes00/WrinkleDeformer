#pragma once

#include <maya/MItMeshEdge.h>
#include <maya/MFnMesh.h>
#include <maya/MStatus.h>
#include <maya/MGlobal.h>
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>



struct Vertex {
	std::vector<uint32_t> vert_faces;
	std::vector<uint32_t> vert_edges;
	std::vector<uint32_t> connected_verts;
};


class meshTopology {
	
	std::vector<Vertex> mesh_verts;
	std::vector<float> rest_edgeLengths;


public:
	meshTopology() {};
	~meshTopology() {};
	void buildFromMesh(MObject& mesh);
	std::map<std::pair<int, int>, float> vertstoRestLen;
	std::unordered_map<int, std::vector<int>> vertIDConnections;

	//const Vertex& getVertex(uint32_t index) const;
};
