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
	
	unsigned int FindAddVertice(const Vector3d& coord, PolyhedralMesh& meshTriangulated, unsigned int& k1, unsigned int k3)
	{
		double tol = 1e-12;

		for (unsigned int i = 0; i < meshTriangulated.Cell0DsCoordinates.cols(); i++) {
			if ((meshTriangulated.Cell0DsCoordinates.col(i) - coord).norm() < tol) {
				meshTriangulated.Cell2DsVertices[k3].push_back(i);
				return i;
			}
		}

		// Nuovo vertice
		meshTriangulated.Cell0DsCoordinates.col(k1) = coord;
		meshTriangulated.Cell0DsId[k1] = k1;
		meshTriangulated.Cell2DsVertices[k3].push_back(k1);
		return k1++;
	}
	
	
	unsigned int FindAddEdge(unsigned int a, unsigned int b, PolyhedralMesh& meshTriangulated, unsigned int& k2, unsigned int k3)
	{
		// Controllo se l'edge (a,b) o (b,a) esiste già
		for (unsigned int i = 0; i < k2; ++i) {
			if ((meshTriangulated.Cell1DsExtrema(i, 0) == a && meshTriangulated.Cell1DsExtrema(i, 1) == b) ||
				(meshTriangulated.Cell1DsExtrema(i, 0) == b && meshTriangulated.Cell1DsExtrema(i, 1) == a)) {

				// Edge già presente ⇒ associarlo alla faccia
				meshTriangulated.Cell2DsEdges[k3].push_back(i);
				return i; // restituisci l'ID dell'edge
			}
		}

		// Edge non presente ⇒ lo aggiungo
		meshTriangulated.Cell1DsExtrema(k2, 0) = a;
		meshTriangulated.Cell1DsExtrema(k2, 1) = b;
		meshTriangulated.Cell1DsId[k2] = k2;

		meshTriangulated.Cell2DsEdges[k3].push_back(k2);
		return k2++; // restituisci il vecchio valore di k2, poi incrementalo
	}
	
	
	/*
	void FindAddVertice(Vector3d coord, PolyhedralMesh& meshTriangulated, unsigned int& k1, unsigned int k3){
		double tol = 1e-12;
		bool found = false;
		
		for (unsigned int i = 0; i< meshTriangulated.Cell0DsCoordinates.cols(); i++){
			
			if ((meshTriangulatedCell0Ds.Coordinates.col(i) - coord).norm() < tol){
				found = true;
				// i vertici esitono già, li aggiungo alla faccia k3 in posizione i
				meshTriangulated.Cell2DsVertices[k3].push_back(i);
				break;
			}
		}
		
		if (!found) {
			// i vertici non esitono, li aggiungo alle coordinate
			meshTriangulated.Cell0DsCoordinates.col(k1)=coord;
			meshTriangulated.Cell0DsId[k1]=k1;
			// aggiungo i vertici alla faccia k3 
			meshTriangulated.Cell2DsVertices[k3].push_back(k1);
			k1++;
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
		
		if (!found) {
			// Aggiungo nuovo lato
			meshTriangulated.Cell1DsExtrema(k2, 0) = a;
			meshTriangulated.Cell1DsExtrema(k2, 1) = b;
			meshTriangulated.Cell1DsId[k2] = k2;

			// Aggiungo alla faccia corrente
			meshTriangulated.Cell2DsEdges[k3].push_back(k2);

			k2++;
		}
	}
	*/	
		
		
		
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
						break;
					}
				}
				if (flag != maxFlag){
					
					Vector3d mediumPoint = (mesh.Cell0DsCoordinates.col(faceVertices[e])+mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]))/2.0;


					// TRIANGOLO SOPRA IL PUNTO MEDIO
					// id della faccia k3 si trova in posizione k3
					meshTriangulated.Cell2DsId[k3]=k3;

					
					// adesso se voglio risalire all'id del punto che ho aggiunto come faccio? 
					// conviene che la funzione restituisca l'id del punto, sia che già esista sia che io lo aggiunga?
					

					unsigned int idA = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, k1, k3);
					unsigned int idB = FindAddVertice(mediumPoint, meshTriangulated, k1, k3);
					unsigned int idC = FindAddVertice(barycenter, meshTriangulated, k1, k3);
					
					// lato vertice - punto medio
					FindAddEdge(idA, idB, meshTriangulated, k2, k3);
					
					// lato punto medio- baricentro
					FindAddEdge(idB, idC, meshTriangulated, k2, k3);
					
					// lato baricentro - vertice
					FindAddEdge(idC, idA, meshTriangulated, k2, k3);
					
					k3++
					
					/*
					FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, k1, k3);
					FindAddVertice(mediumPoint, meshTriangulated, k1, k3);
					FindAddVertice(barycenter, meshTriangulated, k1, k3);
					
					//meshTriangulatedCell2DsVertices[k3]= {k1-2, k1-1, k1};
					
					// lato vertice - punto medio
					FindAddEdge(mesh.Cell0DsCoordinates.col(faceVertices[e]), mediumPoint, meshTriangulated, k2, k3);
					
					// lato punto medio- baricentro
					FindAddEdge(mediumPoint, mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, k2, k3);
					
					// lato baricentro - vertice
					FindAddEdge(barycenter, k1-2, meshTriangulated, k2, k3);
					k3++;
					*/
					
					// triangolo sotto il punto medio
					
					// id della faccia k3 si trova in posizione k3
					meshTriangulated.Cell2DsId[k3]=k3;
					
					unsigned int idA1 = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e+1)%3]), meshTriangulated, k1, k3);
					unsigned int idB1 = FindAddVertice(mediumPoint, meshTriangulated, k1, k3);
					unsigned int idC1 = FindAddVertice(barycenter, meshTriangulated, k1, k3);
					
					// lato vertice - punto medio
					FindAddEdge(idA1, idB1, meshTriangulated, k2, k3);
					
					// lato punto medio- baricentro (questo dovrebbe già esistere)
					FindAddEdge(idB1, idC1, meshTriangulated, k2, k3);
					
					// lato baricentro - vertice
					FindAddEdge(idC1, idA1, meshTriangulated, k2, k3);
					
					k3++
					
					
					/*
					meshTriangulatedCell0DsCoordinates[k1]= mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]);
					meshTriangulatedCell0DsId[k1]= k1;
					k1++;
					
					meshTriangulated.Cell1DsExtrema[k2] << k1, k1-2;
					meshTriangulated.Cell1DsId[k2]= k2;
					k2++;
					meshTriangulated.Cell1DsExtrema[k2] << k1, k1-1;
					meshTriangulated.Cell1DsId[k2]= k2;
					k2++;
					
					meshTriangulated.Cell2DsVertices[k3]= {k1-2, k1-1, k1};
					meshTriangulated.Cell2DsEdges[k3]= {k2-3, k2-1, k2};
					meshTriangulated.Cell2DsId[k3]=k3;
					k3++;
					*/
					
					

				
			} 
			
	
        

        }

