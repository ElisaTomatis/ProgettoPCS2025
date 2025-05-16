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
		int F = dimension[2];
		
		if (q == 3) {
			V += 2*V + E *(subdivisionLevel - 1);
			E += E*subdivisionLevel;
		}
		else if (q == 4) {
			V += 2*V + E *(subdivisionLevel - 1);
			E += E*subdivisionLevel;
		}
		else {
			V += 2*V + E *(subdivisionLevel - 1);
			V += V * subdivisionLevel;
		}
		result[0] = V;  
		result[1] = E;  
		result[2] = F;
		
		return result;
	}
	
	void RemoveDuplicatedVertices(
    Eigen::MatrixXd& Cell0DsCoordinates,
    vector<vector<unsigned int>>& Cell0DsFlag)
	{
		double tol = 1e-12;
		unsigned int maxFlag = numeric_limits<unsigned int>::max();
		size_t n = Cell0DsCoordinates.cols(); // Numero di vertici
	
		for (size_t i = 0; i < n; ++i) {
			if (Cell0DsFlag[i][0] == maxFlag)
				continue;
	
			for (size_t j = i + 1; j < n; ++j) {
				if (Cell0DsFlag[j][0] == maxFlag)
					continue;
	
				// Verifica se condividono almeno un lato
				bool commonSide = false;
				for (unsigned int fi : Cell0DsFlag[i]) {
					for (unsigned int fj : Cell0DsFlag[j]) {
						if (fi == fj) {
							commonSide = true;
							break;
						}
					}
					if (commonSide) break;
				}
	
				if (commonSide) {
					if ((Cell0DsCoordinates.col(i) - Cell0DsCoordinates.col(j)).norm() < tol) {
						Cell0DsFlag[j] = {maxFlag};
					}
				}
			}
		}
	}
	
	void RemoveDuplicatedEdges(
    Eigen::MatrixXi& Cell1DsExtrema,
    vector<vector<unsigned int>>& Cell1DsFlag)
	{
		unsigned int maxFlag = std::numeric_limits<unsigned int>::max();
		size_t n = Cell1DsExtrema.rows(); // Numero di lati
	
		for (size_t i = 0; i < n; ++i) {
			if (Cell1DsFlag[i][0] == maxFlag)
				continue;
	
			for (size_t j = i + 1; j < n; ++j) {
				if (Cell1DsFlag[j][0] == maxFlag)
					continue;
	
				bool commonFlag = false;
				for (unsigned int fi : Cell1DsFlag[i]) {
					for (unsigned int fj : Cell1DsFlag[j]) {
						if (fi == fj) {
							commonFlag = true;
							break;
						}
					}
					if (commonFlag) break;
				}
	
				if (commonFlag){
					int i0 = Cell1DsExtrema(i, 0);
					int i1 = Cell1DsExtrema(i, 1);
					int j0 = Cell1DsExtrema(j, 0);
					int j1 = Cell1DsExtrema(j, 1);
		
					if ((i0 == j0 && i1 == j1) || (i0 == j1 && i1 == j0)){
						Cell1DsFlag[j] = {maxFlag};
					}
			}
		}
	}
}
}