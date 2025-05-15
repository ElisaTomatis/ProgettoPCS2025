#include <vector>
#include <array>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Funzione principale per la triangolazione e salvataggio nella mesh
void triangulateAndStore(PolyhedralMesh& mesh, PolyhedralMesh meshTriangulated, unsigned int b, unsigned int c, const Vector3i& dimension) {
	
	unsigned int subdivisionLevel = b + c;
	
	meshTriangulated.Cell0DsId.reserve(dimension[0]);
	meshTriangulated.Cell0DsCoordinates.reserve(dimension[0],3);
	meshTriangulated.Cell0DsFlag.reserve(dimension[0]);
	
	meshTriangulated.Cell1DsId.reserve(dimension[1]);
	meshTriangulated.Cell1DsExtrema.reserve(dimension[1],2);
	meshTriangulated.Cell1DsFlag.reserve(dimension[1]);
	
	meshTriangulated.Cell2DsId.reserve(dimension[2]);
	meshTriangulated.Cell2DsVertices.reserve(dimension[2],3);
	meshTriangulated.Cell2DsEdges.reserve(dimension[2],3);
	
	unsigned int k1=0;
	unsigned int k2=0;
	unsigned int k3=0;
	
    for (unsigned int faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId) {
        const auto& face = mesh.Cell2DsVertices[faceId];
        Vector2d V0 = mesh.Cell0DsCoordinates.row(face[0]);
        Vector2d V1 = mesh.Cell0DsCoordinates.row(face[1]);
        Vector2d V2 = mesh.Cell0DsCoordinates.row(face[2]);

        vector<vector<int>> vertexGrid;

        for (int i = 0; i <= subdivisionLevel; ++i) {
	        vector<int> row;
            Vector2d start = ((double)i / subdivisionLevel) * V1 + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0;
            Vector2d end   = ((double)i / subdivisionLevel) * V2 + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0;

            for (int j = 0; j <= i; ++j) {
                Vector3d point;
                if (i == 0) {
                    point = V0;
                } else {
                    point = ((double)j / i) * end + ((double)(i - j) / i) * start;
                }

    			meshTriangulated.Cell0DsCoordinates.row(k1) = point.transpose();
    			mesh.Triangulated.Cell0DsId(k1)=k1;
    			
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
					meshTriangulated.Cell0DsFlag[k1] = {-1};
				}

    			k1++;
                row.push_back(k1);
            }
            vertexGrid.push_back(row);
        }

        for (int i = 0; i < subdivisionLevel; ++i) {
            for (int j = 0; j < i; ++j) {
                int v1 = vertexGrid[i][j];
                int v2 = vertexGrid[i + 1][j];
                int v3 = vertexGrid[i + 1][j + 1];

                array<int, 3> verts1 = {v1, v2, v3};
                meshTriangulated.Cell2DsVertices.push_back(verts1);
                meshTriangulated.Cell2DsId(k3)=k3;
                
                for (int e = 0; e < 3; ++e) {
                    int a = verts1[e];
                    int b = verts1[(e + 1) % 3];
                    FindAddEdge(a, b, meshTriangulated, k2, k3)
		        }
		        
		        k3++;
                    
                int v4 = vertexGrid[i][j];
                int v5 = vertexGrid[i + 1][j + 1];
                int v6 = vertexGrid[i][j + 1];

                array<int, 3> verts2 = {v4, v5, v6};
                meshTriangulated.Cell2DsVertices.push_back(verts2);
                meshTriangulated.Cell2DsId(k3)=k3;
                
                for (int e = 0; e < 3; ++e) {
                    int a = verts2[e];
                    int b = verts2[(e + 1) % 3];
                    FindAddEdge(a, b, meshTriangulated, k2, k3)
		            }
                }

            // Triangolo finale in basso a sinistra
            int v1 = vertexGrid[i][i];
            int v2 = vertexGrid[i + 1][i];
            int v3 = vertexGrid[i + 1][i + 1];

            array<int, 3> verts = {v1, v2, v3};
            meshTriangulated.Cell2DsVertices.push_back(verts);
            meshTriangulated.Cell2DsId(k3)=k3;
            
            for (int e = 0; e < 3; ++e) {
                int a = verts[e];
                int b = verts[(e + 1) % 3];
                FindAddEdge(a, b, meshTriangulated, k2, k3)
             }
             k3++;
        }
    }
}


void FindAddEdge(
    unsigned int a, unsigned int b,
    PolyhedralMesh& meshTriangulated,
    unsigned int& k2,
    unsigned int k3)
{
    bool found = false;

	for (int i = 0; i < k2; ++i) {  // Scorri solo fino all'ultimo edge inserito
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
		meshTriangulated.Cell1DsId(k2) = k2;
	
		const auto& c = meshTriangulated.Cell0DsFlag[a];
		const auto& d = meshTriangulated.Cell0DsFlag[b];
	
		bool common = false;
		for (int s = 0; s < c.size(); ++s) {
			for (int t = 0; t < d.size(); ++t) {
				if (c[s] == d[t]) {
					meshTriangulated.Cell1DsFlag[k2] = {c[s]};
					common = true;
					break;
				}
			}
			if (common) break;
		}
	
		if (!common) {
			meshTriangulated.Cell1DsFlag[k2] = {-1};
		}
	
		meshTriangulated.Cell2DsEdges[k3].push_back(k2);
		++k2;
	}
}

