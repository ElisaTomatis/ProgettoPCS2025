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
    void triangulateAndStore2(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated, unsigned int b, unsigned int c, const vector<int>& dimensionDuplicated) {

        unsigned int subdivisionLevel = b+c;
        unsigned int maxFlag = numeric_limits<unsigned int>::max();

        meshTriangulated.Cell0DsId.resize(50);
        meshTriangulated.Cell0DsCoordinates = MatrixXd::Zero(3, 50);
        meshTriangulated.Cell0DsFlag.resize(50);

        meshTriangulated.Cell1DsId.resize(50);
        meshTriangulated.Cell1DsExtrema = MatrixXi::Zero(50, 2);
        meshTriangulated.Cell1DsFlag.resize(50);

        meshTriangulated.Cell2DsId.resize(50);
        meshTriangulated.Cell2DsVertices.resize(50);
        meshTriangulated.Cell2DsEdges.resize(50);

        unsigned int k1 = 0; // Contatore nodi globali (0D)
        unsigned int k2 = 0; // Contatore edge globali (1D)
        unsigned int k3 = 0; // Contatore facce/triangoli globali (2D)

        for (unsigned int faceId = 0; faceId < mesh.Cell2DsId.size(); ++faceId) {
            const auto& face = mesh.Cell2DsVertices[faceId];
            Vector3d V0_base = mesh.Cell0DsCoordinates.col(face[0]);
            Vector3d V1_base = mesh.Cell0DsCoordinates.col(face[1]);
            Vector3d V2_base = mesh.Cell0DsCoordinates.col(face[2]);

            vector<vector<unsigned int>> vertexGrid;

            for (unsigned int i = 0; i <= subdivisionLevel; ++i) {
                vector<unsigned int> row;
                Vector3d start = ((double)i / subdivisionLevel) * V1_base + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0_base;
                Vector3d end = ((double)i / subdivisionLevel) * V2_base + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0_base;
				
				if (i % 2 == 0) {
					for (unsigned int j = 0; j <= i; ++j) {
						Vector3d point;
						if (i == 0) {
							point = V0_base;
						} else {
							point = ((double)j / i) * end + ((double)(i - j) / i) * start;
						}
	
						meshTriangulated.Cell0DsCoordinates.col(k1) = point;
						meshTriangulated.Cell0DsId[k1] = k1;
	
						if (i == 0) {
							meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][2]};
						} else if (i == subdivisionLevel) {
							if (j == 0)
								meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][1]};
							else if (j == subdivisionLevel)
								meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][1], mesh.Cell2DsEdges[faceId][2]};
							else
								meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][1]};
						} else if (j == 0) {
							meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0]};
						} else if (j == i) {
							meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][2]};
						} else {
							meshTriangulated.Cell0DsFlag[k1] = {maxFlag};
						}
						row.push_back(k1);
						k1++;
					}
                vertexGrid.push_back(row); 
            	}
            else {
	            meshTriangulated.Cell0DsCoordinates.col(k1) = start;
				meshTriangulated.Cell0DsId[k1] = k1;
				meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0]};
				row.push_back(k1);
				k1++;
				meshTriangulated.Cell0DsCoordinates.col(k1) = end;
				meshTriangulated.Cell0DsId[k1] = k1;
				meshTriangulated.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][2]};
				row.push_back(k1);
				vertexGrid.push_back(row);
            }
            }
            for (unsigned int i = 0; i < vertexGrid.size(); ++i) {
				cout << "Riga " << i << ": ";
				for (unsigned int j = 0; j < vertexGrid[i].size(); ++j) {
					cout << vertexGrid[i][j] << " ";
				}
				cout << endl;
			}
			
            for (unsigned int i = 0; i < subdivisionLevel; ++i) {
                for (unsigned int j = 0; j <=i; ++j) {
	                
	                if (i%2 == 0){
						unsigned int v1 = vertexGrid[i][j];
						unsigned int v2 = vertexGrid[i + 1][j];
						unsigned int v3 = vertexGrid[i + 1][j + 1];
	
						Vector3d v1_coords = meshTriangulated.Cell0DsCoordinates.col(v1);
						Vector3d v2_coords = meshTriangulated.Cell0DsCoordinates.col(v2);
						Vector3d v3_coords = meshTriangulated.Cell0DsCoordinates.col(v3);
						meshTriangulated.Cell0DsCoordinates.col(k1) = (v1_coords + v2_coords + v3_coords) / 3.0;
						meshTriangulated.Cell0DsId[k1] = k1;
						meshTriangulated.Cell0DsFlag[k1] = {maxFlag}; // il baricentro è sicuramente interno
						
						// aggiungo i lati con congiungo i baricentri con i vertici originali
						vector<unsigned int> verts1 = {v1, v2, v3};
						for (unsigned int e = 0; e < 3; ++e) {
							int a = verts1[e];
							meshTriangulated.Cell1DsExtrema.row(k2) << k1, a; // baricentro - vertice
							meshTriangulated.Cell1DsId[k2] = k2;
							meshTriangulated.Cell1DsFlag[k2] = maxFlag; // si trova sicuramente in centro
							++k2; 
						}

						unsigned int v4 = vertexGrid[i][j];
						unsigned int v5 = vertexGrid[i + 1][j + 1];
						unsigned int v6 = vertexGrid[i][j + 1];
	
						Vector3d v4_coords = meshTriangulated.Cell0DsCoordinates.col(v4);
						Vector3d v5_coords = meshTriangulated.Cell0DsCoordinates.col(v5);
						Vector3d v6_coords = meshTriangulated.Cell0DsCoordinates.col(v6);
						meshTriangulated.Cell0DsCoordinates.col(k1) = (v4_coords + v5_coords + v6_coords) / 3.0;
						meshTriangulated.Cell0DsId[k1] = k1;
						meshTriangulated.Cell0DsFlag[k1] = {maxFlag}; // il baricentro è sicuramente interno
						
						vector<unsigned int> verts2 = {v3, v4, v5};
						for (unsigned int e = 0; e < 3; ++e) {
							int a = verts2[e];
							meshTriangulated.Cell1DsExtrema.row(k2) << k1, a; // baricentro - vertice
							meshTriangulated.Cell1DsId[k2] = k2;
							meshTriangulated.Cell1DsFlag[k2] = maxFlag; // si trova sicuramente in centro
							++k2;
						}
					}

                }
            }
            
            // Creo una mappa che ad ogni lato orginale	 mi associ gli id dei vertici ai lati
        }
    }
}

