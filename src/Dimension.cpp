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
		vector<int> result(3, -1); // inizializza un vettore con 3 valori, tutti -1 
	
		int T = b * b + b * c + c * c;
		int V, E, F;
	
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
		vector<int> result(3, 1);
		int subdivisionLevel = b + c;
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
}