#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	vector<int> ComputePolyhedronVEF(int q, int b, int c)
	{
		vector<int> result(3); // inizializza un vettore con 3 valori, tutti -1 
	
		int T = 0;
		T = b * b + b * c + c * c;
		int V = 0;
		int E = 0;
		int F = 0;
	
		if (q == 3) {
			V = 2 * T + 2;
			E = 6 * T;
			F = 4 * T;
		}
		else if (q == 4) {
			V = 4 * T + 2;
			E = 12 * T;
			F = 8 * T;
		}
		else {
			V = 10 * T + 2;
			E = 30 * T;
			F = 20 * T;
		}
	
		result[0] = V;  // Primo elemento: V (vertici)
		result[1] = E;  // Secondo elemento: E (spigoli)
		result[2] = F;  // Terzo elemento: F (facce)
	
		return result;  // Restituisce il vettore con i valori di V, E, F
	}
	
	vector<int> CalculateDuplicated(int q, int b, int c, const vector<int>& dimension)
	{
		vector<int> result(3);
		int subdivisionLevel = 0;
		subdivisionLevel = b + c;
		int V = dimension[0];
		int E = dimension[1];
		
		if (q == 3) {
			V += 2*4 + 6*(subdivisionLevel - 1);
			E += 6*subdivisionLevel;
		}
		else if (q == 4) {
			V += 3*6 + 12*(subdivisionLevel - 1);
			E += 12*subdivisionLevel;
		}
		else {
			V += 4*12 + 30*(subdivisionLevel - 1);
			E += 30* subdivisionLevel;
		}
		result[0] = V;  
		result[1] = E;  
		result[2] = dimension[2]; //il numero delle facce non cambia
		
		return result;
	}
	
	void RemoveDuplicatedVertices(PolyhedralMesh& meshTriangulated)
	{
		double tol = 1e-12;
		unsigned int maxFlag = numeric_limits<unsigned int>::max();
		size_t n = meshTriangulated.Cell0DsCoordinates.cols();
		vector<bool> reached(n, true);

		for (size_t i = 0; i < n; ++i) {
			if (meshTriangulated.Cell0DsFlag[i][0] == maxFlag)
				continue;
			
			for (size_t j = i+1; j < n; ++j) {
				if (meshTriangulated.Cell0DsFlag[j][0] == maxFlag)
					continue;
				
				bool commonSide = false;
				for (unsigned int fi : meshTriangulated.Cell0DsFlag[i]) {
					for (unsigned int fj : meshTriangulated.Cell0DsFlag[j]) {
						if (fi == fj) {
							commonSide = true;
							break;
						}
					}
					if (commonSide) break;
				}
					
				if (commonSide) {
					if ((meshTriangulated.Cell0DsCoordinates.col(i) - meshTriangulated.Cell0DsCoordinates.col(j)).norm() < tol) {
						reached[i]=false;
						break;
					}
				}
			}
			if (reached[i]) {
				meshTriangulated.Cell0DsFlag[i] = {maxFlag};
			}
		}
	}
	
	void RemoveDuplicatedEdges(PolyhedralMesh& meshTriangulated)
	{
		double tol = 1e-12;
		unsigned int maxFlag = numeric_limits<unsigned int>::max();
		size_t n = meshTriangulated.Cell1DsExtrema.rows();
		vector<bool> reached(n, false);
			
		for (size_t i = 0; i < n; ++i) {
			if (meshTriangulated.Cell1DsFlag[i] == maxFlag || reached[i]){
				continue; 
			}
			
			for (size_t j = i + 1; j < n; ++j) {
				if (meshTriangulated.Cell1DsFlag[j] == maxFlag || reached[j]){
					continue;
				}
	
				if (meshTriangulated.Cell1DsFlag[i] == meshTriangulated.Cell1DsFlag[j]) {
					int i0 = meshTriangulated.Cell1DsExtrema(i, 0);
					int i1 = meshTriangulated.Cell1DsExtrema(i, 1);
					int j0 = meshTriangulated.Cell1DsExtrema(j, 0);
					int j1 = meshTriangulated.Cell1DsExtrema(j, 1);
					
					if (((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol && (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol) || ((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol && (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol)) {
						meshTriangulated.Cell1DsFlag[i] = maxFlag;
						reached[j]=true;
						break;
					}
				}
			}
		}	
	}
	
	void NewMesh(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension)
	{
		unsigned int maxFlag = std::numeric_limits<unsigned int>::max();
		
		meshFinal.Cell0DsId.resize(dimension[0]);
		meshFinal.Cell0DsCoordinates = MatrixXd::Zero(3, dimension[0]);
		
		meshFinal.Cell1DsId.resize(dimension[1]);
		meshFinal.Cell1DsExtrema = MatrixXi::Zero(dimension[1], 2);
		
		meshFinal.Cell2DsId.resize(dimension[2]);
		meshFinal.Cell2DsVertices.resize(dimension[2]);
		meshFinal.Cell2DsEdges.resize(dimension[2]);
		
		unsigned int k1=0;
		unsigned int k2=0;
		
		// VERTICI
		
		/*
		for (size_t i = 0; i < meshTriangulated.Cell0DsFlag.size(); ++i) {
			for (size_t j = 0; j < meshTriangulated.Cell0DsFlag[i].size(); ++j) {
				if (meshTriangulated.Cell0DsFlag[i][j] == maxFlag){
					meshFinal.Cell0DsCoordinates.col(k1)= meshTriangulated.Cell0DsCoordinates.col(i);
					k1 ++;
				}
			}
		}*/
		
		for (size_t i = 0; i < meshTriangulated.Cell0DsFlag.size(); ++i) {
			bool hasMaxFlag = false;
			for (size_t j = 0; j < meshTriangulated.Cell0DsFlag[i].size(); ++j) {
				if (meshTriangulated.Cell0DsFlag[i][j] == maxFlag) {
					hasMaxFlag = true;
					break; // appena trovato, esci
				}
			}
			if (hasMaxFlag) {
				meshFinal.Cell0DsCoordinates.col(k1) = meshTriangulated.Cell0DsCoordinates.col(i);
				++k1;
			}
		}

		
		for (int i=0; i<dimension[0]; ++i){
			meshFinal.Cell0DsId[i]=i;
		}
		
		// LATI
		for (size_t i = 0; i < meshTriangulated.Cell1DsFlag.size(); ++i) {
			if (meshTriangulated.Cell1DsFlag[i] == maxFlag){
				meshFinal.Cell1DsExtrema.row(k2)= meshTriangulated.Cell1DsExtrema.row(i);
				k2++;
			}	
		}
		
		for (int i=0; i<dimension[1]; ++i){
			meshFinal.Cell1DsId[i]=i;
		}
		
		// FACCE
		for (int i=0; i<dimension[2]; ++i){
			meshFinal.Cell2DsId[i]=i;
		}
		
		meshFinal.Cell2DsVertices = meshTriangulated.Cell2DsVertices;
		meshFinal.Cell2DsEdges = meshTriangulated.Cell2DsEdges;	
	}
}
