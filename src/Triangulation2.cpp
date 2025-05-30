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
	void FindAddVertice(Vector3d coord, PolyhedralMesh& meshTriangulated, unsigned int k1){
		double tol = 1e-12;
		bool found = false;
		
		for (unsigned int i = 0; i< meshTriangulatedCell0DsCoordinates.cols(); i++){
			if ((meshTriangulatedCell0DsCoordinates.col(i) - coord).norm() < tol){
				found = true;
				meshTriangulatedCell2DsVertices[k3].push_back(i);
				break;
			}
		}
				meshTriangulatedCell0DsCoordinates[k1]=coord;
				meshTriangulatedCell0DsId[k1]=k1;
				k1++;
				break;
			}
		}
	}
	
	
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
		
		
		
		
		
    void triangulateAndStore2(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated, unsigned int b, unsigned int c, const vector<int>& dimension) {

        unsigned int subdivisionLevel = b+c;
        unsigned int maxFlag = numeric_limits<unsigned int>::max();

        meshTriangulated.Cell0DsId.resize(dimension[0]);
        meshTriangulated.Cell0DsCoordinates = MatrixXd::Zero(3, dimension[0]);
        meshTriangulated.Cell0DsFlag.resize(dimension[0]);

        meshTriangulated.Cell1DsId.resize(dimension[1]);
        meshTriangulated.Cell1DsExtrema = MatrixXi::Zero(dimension[1], 2);
        meshTriangulated.Cell1DsFlag.resize(dimension[1]);

        meshTriangulated.Cell2DsId.resize(dimension[2]);
        meshTriangulated.Cell2DsVertices.resize(dimension[2]);
        meshTriangulated.Cell2DsEdges.resize(dimension[2]);

        unsigned int k1 = 0; // Contatore nodi globali (0D)
        unsigned int k2 = 0; // Contatore edge globali (1D)
        unsigned int k3 = 0; // Contatore facce/triangoli globali (2D)
        
        for (unsigned int faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId){
			const auto& faceVertices = mesh.Cell2DsVertices[faceId];
			Vector3d V0 = mesh.Cell0DsCoordinates.col(faceVertices[0]);
			Vector3d V1 = mesh.Cell0DsCoordinates.col(faceVertices[1]);
			Vector3d V2 = mesh.Cell0DsCoordinates.col(faceVertices[2]);
			
			Vector3d barycenter = ((V0 + V1 + V2)/3.0);
			
			const auto& faceEdges = mesh.Cell2DsEdges[faceId];
			for (unsigned int e = 0; e <3; e++){
				unsigned int flag = 0;
				
				for(unsigned int i=0; i < mesh.Cell1DsId.size(); i++){
					if (faceEdges[e] == mesh.Cell1DsId[i]){
						flag = mesh.Cell1DsFlag[i];
						break
					}
				}
				if (flag != maxFlag){
					
					// triangolo sopra il punto medio
					Vector3d mediumPoint = (mesh.Cell0DsCoordinates.col(faceVertices[e])+mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]))/2.0;

					FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, k1);
					FindAddVertice(mediumPoint, meshTriangulated, k1);
					FindAddVertice(barycenter, meshTriangulated, k1);
					
					meshTriangulatedCell2DsVertices[k3]= {k1-2, k1-1, k1};
					meshTriangulatedCell2DsId[k3]=k3;
					FindAddEdge(k1-2, k1-1, meshTriangulated, k2, k3);
					FindAddEdge(k1-1, k1, meshTriangulated, k2, k3);
					FindAddEdge(k1, k1-2, meshTriangulated, k2, k3);
					k3++;
					
					// triangolo sotto il punto medio
					meshTriangulatedCell0DsCoordinates[k1]= mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]);
					meshTriangulatedCell0DsId[k1]= k1;
					k1++;
					
					meshTriangulatedCell1DsExtrema[k2] << k1, k1-2;
					meshTriangulatedCell1DsId[k2]= k2;
					k2++;
					meshTriangulatedCell1DsExtrema[k2] << k1, k1-1;
					meshTriangulatedCell1DsId[k2]= k2;
					k2++;
					
					meshTriangulatedCell2DsVertices[k3]= {k1-2, k1-1, k1};
					meshTriangulatedCell2DsEdges[k3]= {k2-3, k2-1, k2};
					meshTriangulatedCell2DsId[k3]=k3;
					k3++;
					
				
			} 
			
	
        

        }

