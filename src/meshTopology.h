#pragma once

#include <maya/MItMeshEdge.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MPointArray.h>
#include <maya/MMatrix.h>
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

struct TriangleData {

	float qInv[2][2];
	int vertIdx[3];
	MVector normal[2];
	int windingSign;
};


class meshTopology {
	
	std::vector<Vertex> mesh_verts;
	std::vector<float> rest_edgeLengths;


public:
	meshTopology() {};
	~meshTopology() {};
	void buildFromMesh(MObject& mesh, int numVerts);
	std::vector<float> vertstoRestLen;
	std::map<std::pair<int, int>, float> map_vertsToRL;
	//std::unordered_map<int, std::vector<int>> vertIDConnections;
	std::vector<int> adjacencyData;    // flat list of all neighbor indices
	std::vector<int> adjacencyStart;   // adjacencyStart[v] = start index in adjacencyData
	std::vector<int> adjacencyCount;   // number of neighbors for vertex v
	std::vector<TriangleData> tritoQInv;

	//const Vertex& getVertex(uint32_t index) const;
};
