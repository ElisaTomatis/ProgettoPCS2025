#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	void FindAddEdge(
		int a, int b,
		PolyhedralMesh& meshTriangulated,
		unsigned int& k2,
		unsigned int k3)
	{
		bool found = false;
	
		for (unsigned int i = 0; i <= k2; ++i) {  // Scorri solo fino all'ultimo edge inserito
			if ((meshTriangulated.Cell1DsExtrema(i, 0) == a && meshTriangulated.Cell1DsExtrema(i, 1) == b) ||
				(meshTriangulated.Cell1DsExtrema(i, 0) == b && meshTriangulated.Cell1DsExtrema(i, 1) == a)) {
				
				// Edge già presente ⇒ usalo
				meshTriangulated.Cell2DsEdges[k3].push_back(i);
				found = true;
				break;
			}
		}
		
		if (!found) {
			// Edge non esiste ⇒ lo creiamo
			meshTriangulated.Cell1DsExtrema.row(k2) << a, b;
			meshTriangulated.Cell1DsId[k2] = k2;
		
			const auto& c = meshTriangulated.Cell0DsFlag[a];
			const auto& d = meshTriangulated.Cell0DsFlag[b];
		
			bool common = false;
			for (size_t s = 0; s < c.size(); ++s) {
				for (size_t t = 0; t < d.size(); ++t) {
					if (c[s] == d[t]) {
						meshTriangulated.Cell1DsFlag[k2] = c[s];
						common = true;
						break;
					}
				}
				if (common) break;
			}
			
			if (!common) {
				meshTriangulated.Cell1DsFlag[k2] = numeric_limits<unsigned int>::max();
			}
		
			meshTriangulated.Cell2DsEdges[k3].push_back(k2);
			++k2;
		}
	}
	
	// Funzione principale per la triangolazione e salvataggio nella mesh
	void triangulateAndStore(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated, unsigned int b, unsigned int c, const vector<int>& dimensionDuplicated) {
		
		unsigned int subdivisionLevel = b + c;
		
		meshTriangulated.Cell0DsId.resize(dimensionDuplicated[0]);
		meshTriangulated.Cell0DsCoordinates = MatrixXd::Zero(3, dimensionDuplicated[0]);
		meshTriangulated.Cell0DsFlag.resize(dimensionDuplicated[0]);
		
		meshTriangulated.Cell1DsId.resize(dimensionDuplicated[1]);
		meshTriangulated.Cell1DsExtrema = MatrixXi::Zero(dimensionDuplicated[1], 2);
		meshTriangulated.Cell1DsFlag.resize(dimensionDuplicated[1]);
		
		meshTriangulated.Cell2DsId.resize(dimensionDuplicated[2]);
		meshTriangulated.Cell2DsVertices.resize(dimensionDuplicated[2]);
		meshTriangulated.Cell2DsEdges.resize(dimensionDuplicated[2]);

		
		/*
		meshTriangulated.Cell2DsVertices.reserve(dimensionDuplicated[2]);
		meshTriangulated.Cell2DsEdges.reserve(dimensionDuplicated[2]);*/
		
		
		unsigned int k1=0;
		unsigned int k2=0;
		unsigned int k3=0;
		
		for (unsigned int faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId){
			const auto& face = mesh.Cell2DsVertices[faceId];
			Vector3d V0 = mesh.Cell0DsCoordinates.col(face[0]);
			Vector3d V1 = mesh.Cell0DsCoordinates.col(face[1]);
			Vector3d V2 = mesh.Cell0DsCoordinates.col(face[2]);
	
			vector<vector<int>> vertexGrid;
	
			for (unsigned int i = 0; i <= subdivisionLevel; ++i) {
				vector<int> row;
				Vector3d start = ((double)i / subdivisionLevel) * V1 + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0;
				Vector3d end   = ((double)i / subdivisionLevel) * V2 + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0;
	
				for (unsigned int j = 0; j <= i; ++j) {
					Vector3d point;
					if (i == 0) {
						point = V0;
					} else {
						point = ((double)j / i) * end + ((double)(i - j) / i) * start;
					}
					
					meshTriangulated.Cell0DsCoordinates.col(k1) = point;
					meshTriangulated.Cell0DsId[k1]=k1;

					if (i == 0) {
					meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][2]};
					}
					else if (i == subdivisionLevel) {
						if (j == 0)
							meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][1]};
						else if (j == subdivisionLevel)
							meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][1], mesh.Cell2DsEdges[faceId][2]};
						else
							meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][1]};
					}
					else if (j == 0) {
						meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0]};
					}
					else if (j == i) {
						meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][2]};
					}
					else {
						meshTriangulated.Cell0DsFlag[k1] = {numeric_limits<unsigned int>::max()};
					}		
	
					row.push_back(k1);
					k1++;
				}
				vertexGrid.push_back(row);
			}
			for (unsigned int i = 0; i < subdivisionLevel; ++i) {
				for (unsigned int j = 0; j < i; ++j) {
					unsigned int v1 = vertexGrid[i][j];
					unsigned int v2 = vertexGrid[i + 1][j];
					unsigned int v3 = vertexGrid[i + 1][j + 1];
		
					vector<unsigned int> verts1 = {v1, v2, v3};
						
					// meshTriangulated.Cell2DsVertices.push_back(verts1);
					meshTriangulated.Cell2DsVertices[k3]=verts1;
						
					meshTriangulated.Cell2DsId[k3]=k3;
						
					for (unsigned int e = 0; e < 3; ++e) {
						int a = verts1[e];
						int b = verts1[(e + 1) % 3];
						FindAddEdge(a, b, meshTriangulated, k2, k3);
					}
	
					k3++;
							
					unsigned int v4 = vertexGrid[i][j];
					unsigned int v5 = vertexGrid[i + 1][j + 1];
					unsigned int v6 = vertexGrid[i][j + 1];
		
					vector<unsigned int> verts2 = {v4, v5, v6};
					//meshTriangulated.Cell2DsVertices.push_back(verts2);
					meshTriangulated.Cell2DsVertices[k3]=verts2;
					meshTriangulated.Cell2DsId[k3]=k3;
						
					for (unsigned int e = 0; e < 3; ++e) {
						int a = verts2[e];
						int b = verts2[(e + 1) % 3];
						FindAddEdge(a, b, meshTriangulated, k2, k3);
					}
					k3++;
				}
		
				// Triangolo finale in basso a sinistra
				unsigned int v1 = vertexGrid[i][i];
				unsigned int v2 = vertexGrid[i + 1][i];
				unsigned int v3 = vertexGrid[i + 1][i + 1];
	
				vector<unsigned int> verts = {v1, v2, v3};
				meshTriangulated.Cell2DsVertices[k3]=verts;
				// meshTriangulated.Cell2DsVertices.push_back(verts);
				meshTriangulated.Cell2DsId[k3]=k3;
					
				for (unsigned int e = 0; e < 3; ++e) {
					int a = verts[e];
					int b = verts[(e + 1) % 3];
					FindAddEdge(a, b, meshTriangulated, k2, k3);
				}
				k3++;
		}
			
		}

	}
	
	void PopulateCell3D(PolyhedralMesh& meshTriangulated, const vector<int>& dimension){
		
		unsigned int maxFlag = numeric_limits<unsigned int>::max();
		meshTriangulated.Cell3DsId = {0};
		meshTriangulated.NumCells0Ds = dimension[0];
		meshTriangulated.NumCells1Ds = dimension[1];
		meshTriangulated.NumCells2Ds = dimension[2];

		for (unsigned int i = 0; i < meshTriangulated.Cell0DsId.size(); i++){
			if (meshTriangulated.Cell0DsFlag[i][0] == maxFlag){
				meshTriangulated.Cell3DsVertices.push_back(meshTriangulated.Cell0DsId[i]);
			}
		}
		for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++){
			if (meshTriangulated.Cell1DsFlag[i] == maxFlag){
				meshTriangulated.Cell3DsEdges.push_back(meshTriangulated.Cell1DsId[i]);
			}
		}
		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); i++){
			meshTriangulated.Cell3DsFaces.push_back(i);
		}
	}
	

}
