#include <meshTopology.h>


void meshTopology::buildFromMesh(MObject& mesh, int numVerts) {


	MStatus status;
	MItMeshEdge edgeIter(mesh, &status);
	MItMeshPolygon polyIter(mesh, &status);
	adjacencyCount.assign(numVerts, 0);

	for (; !edgeIter.isDone(); edgeIter.next()) {
		// get vertex indices for the edge
		int vertId1 = edgeIter.index(0);
		int vertId2 = edgeIter.index(1);

		// store connectivity for each vertex
		adjacencyCount[vertId1]++;
		adjacencyCount[vertId2]++;

		MPoint pointId1 = edgeIter.point(0, MSpace::kObject, &status);
		MPoint pointId2 = edgeIter.point(1, MSpace::kObject, &status);

		double dist = pointId1.distanceTo(pointId2);

		std::pair pointPair = std::make_pair(std::min(vertId1, vertId2), std::max(vertId1, vertId2));
		map_vertsToRL[pointPair] = dist;
		vertstoRestLen.push_back(dist);
	}
	adjacencyStart.resize(numVerts + 1, 0);
	for (int v = 0; v < numVerts; ++v) {
		adjacencyStart[v + 1] = adjacencyStart[v] + adjacencyCount[v];
	}

	// Total neighbors is the last entry
	int totalNeighbors = adjacencyStart[numVerts];
	adjacencyData.resize(totalNeighbors);

	std::vector<int> writeCursor = adjacencyStart; // copy start positions

	for (edgeIter.reset(); !edgeIter.isDone(); edgeIter.next()) {
		int v1 = edgeIter.index(0);
		int v2 = edgeIter.index(1);

		adjacencyData[writeCursor[v1]++] = v2;
		adjacencyData[writeCursor[v2]++] = v1;
	}

	MFnMesh restMesh(mesh);

	// iterate over all the faces in the mesh
	for (; !polyIter.isDone(); polyIter.next()) {
		MPointArray trianglePoints;
		MIntArray triangleVertexIndices;

		status = polyIter.getTriangles(trianglePoints, triangleVertexIndices, MSpace::kObject);

		if (MS::kSuccess == status) {

			for (unsigned int i = 0; i < triangleVertexIndices.length(); i += 3) {
				int v0_idx = triangleVertexIndices[i];
				int v1_idx = triangleVertexIndices[i + 1];
				int v2_idx = triangleVertexIndices[i + 2];

				// get world rest length positions
				MPoint p0 = trianglePoints[i];
				MPoint p1 = trianglePoints[i + 1];
				MPoint p2 = trianglePoints[i + 2];

				float2 uv0, uv1, uv2;
				restMesh.getUVAtPoint(p0, uv0, MSpace::kObject);
				restMesh.getUVAtPoint(p1, uv1, MSpace::kObject);
				restMesh.getUVAtPoint(p2, uv2, MSpace::kObject);

				MVector uv_0(uv0[0], uv0[1]);
				MVector uv_1(uv1[0], uv1[1]);
				MVector uv_2(uv2[0], uv2[1]);

				// compute uv space tangents
				// (tu,tv) = (x1 −x0,x2 −x0) (u1 −u0,u2 −u0)^-1
				MVector e1_world = (p1 - p0);  // 3x2
				MVector e2_world = (p2 - p0);
				MVector e1_uv = (uv_1 - uv_0); // 2x2
				MVector e2_uv = (uv_2 - uv_0);

				// [[a, b], [c, d]], inverse = 1/(ad-bc) * [[d, -b], [-c, a]]
				float a = e1_uv[0];
				float b = e2_uv[0];
				float c = e1_uv[1];
				float d = e2_uv[1];
				float det = a * d - b * c;
				float inv_det = 1.0f / det;

				// gotcha: MVector * MVector is dot product, not element-wise multiplication, need to multiply by scalars manually
				MVector tu = e1_world * (d * inv_det) + e2_world * (-c * inv_det);
				MVector tv = e1_world * (-b * inv_det) + e2_world * (a * inv_det);

				MVector n1 = tu.normal();
				MVector ortho_tv = (tu ^ (tu ^ tv)).normal();
				MVector n2 = ortho_tv.normal();

				// take the transpose
				float q00 = n1 * e1_world;
				float q01 = n1 * e2_world;
				float q10 = n2 * e1_world;
				float q11 = n2 * e2_world;

				det = q00 * q11 - q01 * q10;
				inv_det = 1.0f / det;

				// Q⁻¹ = (1/det) * [[ q11, -q01], [-q10,  q00]]
				float qi00 = q11 * inv_det;
				float qi01 = -q01 * inv_det;
				float qi10 = -q10 * inv_det;
				float qi11 = q00 * inv_det;

				TriangleData triData;
				triData.windingSign = (det < 0) ? -1.0f : 1.0f;
				triData.qInv[0][0] = qi00;
				triData.qInv[0][1] = qi01;
				triData.qInv[1][0] = qi10;
				triData.qInv[1][1] = qi11;
				triData.vertIdx[0] = v0_idx;
				triData.vertIdx[1] = v1_idx;
				triData.vertIdx[2] = v2_idx;
				triData.normal[0] = n1;
				triData.normal[1] = n2;
		
				tritoQInv.push_back(triData);
			}
		}
	}
			

	MGlobal::displayInfo(MString("Stored ") + vertstoRestLen.size() + MString(" edges"));


}