#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
std::vector<int> ComputePolyhedronVEF(unsigned int q, int b, int c)
{
    std::vector<int> result(3, -1); // inizializza un vettore con 3 valori, tutti -1 (indicando errore)

    if (q != 3 && q != 4 && q != 5) {
        return result; // Se q non è valido, restituisce il vettore con errori (-1, -1, -1)
    }

    int T = b * b + b * c + c * c;
    unsigned int V, E, F;

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

bool generateTetrahedron(PolyhedralMesh& mesh) {
	
	// VERTICI
    double r = sqrt(3.0) / 3.0;
    
    mesh.Cell0DsCoordinates.reserve(4,3);
    mesh.Cell0DsId.reserve(4);

    mesh.Cell0DsCoordinates(0, 0) = r;   mesh.Cell0DsCoordinates(1, 0) = r;   mesh.Cell0DsCoordinates(2, 0) = r;  
    mesh.Cell0DsCoordinates(0, 1) = -r;  mesh.Cell0DsCoordinates(1, 1) = -r;  mesh.Cell0DsCoordinates(2, 1) = r; 
    mesh.Cell0DsCoordinates(0, 2) = -r;  mesh.Cell0DsCoordinates(1, 2) = r;   mesh.Cell0DsCoordinates(2, 2) = -r;
    mesh.Cell0DsCoordinates(0, 3) = r;   mesh.Cell0DsCoordinates(1, 3) = -r;  mesh.Cell0DsCoordinates(2, 3) = -r;
    
    mesh.Cell0DsId = {0,1,2,3};
    
    // LATI
    mesh.Cell1DsId.reserve(6);
    Cell1DsId = {0,1,2,3,4,5};

    mesh.Cell1DsExtrema.resize(6, 2);
    mesh.Cell1DsExtrema <<
        0, 1,
        1, 2,
        2, 0,
        0, 3,
        3, 1,
        2, 3;

    // FACCE
    mesh.Cell2DsId.reserve(4);
    mesh.Cell2DsId = {0, 1, 2, 3};
    
    mesh.Cell2DsVertices.reserve(4,3);
    mesh.Cell2DsEdges.reserve(4,3);

    mesh.Cell2DsVertices = {
        {0, 1, 2}, 
        {0, 3, 1}, 
        {1, 3, 2},
        {2, 3, 0} 
    };

    mesh.Cell2DsEdges = {
        {0, 1, 2}, 
        {3, 4, 0},
        {4, 5, 1},
        {5, 3, 2} 
    };

    // POLIEDRO
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(4);
    mesh.Cell3DsEdges.reserve(6);
    mesh.Cell3DsFaces.reserve(4);

    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3};
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5};
    mesh.Cell3DsFaces = {0, 1, 2, 3};

    return true;
}

bool generateCube(PolyhedralMesh& mesh) {
	
	// VERTICI
	double r = sqrt(3.0) / 3.0;

	mesh.Cell0DsCoordinates.reserve(8, 3); 
	mesh.Cell0DsId.reserve(8);

	mesh.Cell0DsCoordinates(0, 0) = r;   mesh.Cell0DsCoordinates(1, 0) = r;   mesh.Cell0DsCoordinates(2, 0) = r;  
	mesh.Cell0DsCoordinates(0, 1) = r;   mesh.Cell0DsCoordinates(1, 1) = r;   mesh.Cell0DsCoordinates(2, 1) = -r;  
	mesh.Cell0DsCoordinates(0, 2) = r;   mesh.Cell0DsCoordinates(1, 2) = -r;  mesh.Cell0DsCoordinates(2, 2) = r;  
	mesh.Cell0DsCoordinates(0, 3) = r;   mesh.Cell0DsCoordinates(1, 3) = -r;  mesh.Cell0DsCoordinates(2, 3) = -r;  
	mesh.Cell0DsCoordinates(0, 4) = -r;  mesh.Cell0DsCoordinates(1, 4) = r;   mesh.Cell0DsCoordinates(2, 4) = r;  
	mesh.Cell0DsCoordinates(0, 5) = -r;  mesh.Cell0DsCoordinates(1, 5) = r;   mesh.Cell0DsCoordinates(2, 5) = -r;  
	mesh.Cell0DsCoordinates(0, 6) = -r;  mesh.Cell0DsCoordinates(1, 6) = -r;  mesh.Cell0DsCoordinates(2, 6) = r;  
	mesh.Cell0DsCoordinates(0, 7) = -r;  mesh.Cell0DsCoordinates(1, 7) = -r;  mesh.Cell0DsCoordinates(2, 7) = -r;
	
	mesh.Cell0DsId = {0,1,2,3,4,5,6,7};
	
	// LATI
    mesh.Cell1DsId.reserve(12);
    Cell1DsId = {0,1,2,3,4,5,6,7,8,9,10,11};

    mesh.Cell1DsExtrema.resize(12, 2);
	mesh.Cell1DsExtrema <<
		0, 1,  // Lato 0
		1, 2,  // Lato 1
		2, 3,  // Lato 2
		3, 0,  // Lato 3
		4, 5,  // Lato 4
		5, 6,  // Lato 5
		6, 7,  // Lato 6
		7, 4,  // Lato 7
		0, 4,  // Lato 8
		1, 5,  // Lato 9
		2, 6,  // Lato 10
		3, 7;  // Lato 11
		
	// FACCE
	mesh.Cell2DsId.reserve(6); 
	mesh.Cell2DsId = {0, 1, 2, 3, 4, 5};
	
	mesh.Cell2DsVertices.reserve(6, 4);
	mesh.Cell2DsEdges.reserve(6, 4); 
	
	mesh.Cell2DsVertices = {
		{0, 1, 2, 3},  
		{4, 5, 6, 7},  
		{0, 1, 5, 4},  
		{1, 2, 6, 5},  
		{2, 3, 7, 6},
		{3, 0, 4, 7} 
    };
    
	mesh.Cell2DsEdges = {
		{0, 1, 2, 3},
		{4, 5, 6, 7},
		{0, 9, 4, 8},
		{1, 10, 5, 9},
		{2, 11, 6, 10},
		{3, 8, 7, 11}
	};
	
	// POLIEDRO 
	mesh.Cell3DsId.reserve(1);
	mesh.Cell3DsVertices.reserve(8); 
	mesh.Cell3DsEdges.reserve(12);
	mesh.Cell3DsFaces.reserve(6); 
	
	mesh.Cell3DsId = {0};
	mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5, 6, 7};
	mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5};
	
	return true;
}
	
		
bool generateOctahedron(PolyhedralMesh& mesh) {
	
    // VERTICI
    double r = 1.0; 

    mesh.Cell0DsCoordinates.reserve(6, 3);
    mesh.Cell0DsId.reserve(6);

    mesh.Cell0DsCoordinates(0, 0) = r;   mesh.Cell0DsCoordinates(1, 0) = 0.0;  mesh.Cell0DsCoordinates(2, 0) = 0.0;  
    mesh.Cell0DsCoordinates(0, 1) = -r;  mesh.Cell0DsCoordinates(1, 1) = 0.0;  mesh.Cell0DsCoordinates(2, 1) = 0.0;  
    mesh.Cell0DsCoordinates(0, 2) = 0.0;  mesh.Cell0DsCoordinates(1, 2) = r;   mesh.Cell0DsCoordinates(2, 2) = 0.0;  
    mesh.Cell0DsCoordinates(0, 3) = 0.0;  mesh.Cell0DsCoordinates(1, 3) = -r;  mesh.Cell0DsCoordinates(2, 3) = 0.0;  
    mesh.Cell0DsCoordinates(0, 4) = 0.0;  mesh.Cell0DsCoordinates(1, 4) = 0.0;  mesh.Cell0DsCoordinates(2, 4) = r;   
    mesh.Cell0DsCoordinates(0, 5) = 0.0;  mesh.Cell0DsCoordinates(1, 5) = 0.0;  mesh.Cell0DsCoordinates(2, 5) = -r;  

    mesh.Cell0DsId = {0, 1, 2, 3, 4, 5};
    
    // LATI
    mesh.Cell1DsId.reserve(12);
    Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    mesh.Cell1DsExtrema.resize(12, 2);
    mesh.Cell1DsExtrema <<
        0, 1,  // Lato 0
        0, 3,  // Lato 1
        0, 2,  // Lato 2
        0, 4,  // Lato 3
        5, 4,  // Lato 4
        2, 5,  // Lato 5
        3, 5,  // Lato 6
        1, 5,  // Lato 7
        1, 3,  // Lato 8
        1, 2,  // Lato 9
        4, 2,  // Lato 10
        4, 3;  // Lato 11

    // FACCE
    mesh.Cell2DsId.reserve(8);
    mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7};

    mesh.Cell2DsVertices.reserve(8, 3);
    mesh.Cell2DsEdges.reserve(8, 3);

    mesh.Cell2DsVertices = {
        {0, 1, 3},  // Faccia 0
        {0, 1, 2},  // Faccia 1
        {0, 2, 4},  // Faccia 2
        {0, 3, 4},  // Faccia 3
        {2, 5, 4},  // Faccia 4
        {2, 5, 1},  // Faccia 5
        {1, 3, 5},  // Faccia 6
        {3, 5, 4}   // Faccia 7
    };

    mesh.Cell2DsEdges = {
        {0, 8, 1},  // Faccia 0
        {0, 9, 2},  // Faccia 1
        {2, 10, 3}, // Faccia 2
        {1, 11, 3},  // Faccia 3
        {5, 4, 10}, // Faccia 4
        {5, 7, 9}, // Faccia 5
        {8, 6, 7},  // Faccia 6
        {6, 4, 11}  // Faccia 7
    };

    // POLIEDRO
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(6);
    mesh.Cell3DsEdges.reserve(12);
    mesh.Cell3DsFaces.reserve(8);

    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5};
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7};
    
    return true;
}

bool generateDodecahedron(PolyhedralMesh& mesh) {
	
    // VERTICI
    double r = 1.0;
    double phi = (1.0 + sqrt(5.0)) / 2.0;

    mesh.Cell0DsCoordinates.reserve(20, 3);
    mesh.Cell0DsId.reserve(20);
    double norm = sqrt(3.0);
    
    mesh.Cell0DsCoordinates(0, 0) = r / norm;  mesh.Cell0DsCoordinates(1, 0) = r / norm;  mesh.Cell0DsCoordinates(2, 0) = r / norm;
	mesh.Cell0DsCoordinates(0, 1) = r / norm;  mesh.Cell0DsCoordinates(1, 1) = r / norm;  mesh.Cell0DsCoordinates(2, 1) = -r / norm;
	mesh.Cell0DsCoordinates(0, 2) = r / norm;  mesh.Cell0DsCoordinates(1, 2) = -r / norm;  mesh.Cell0DsCoordinates(2, 2) = r / norm;
	mesh.Cell0DsCoordinates(0, 3) = r / norm;  mesh.Cell0DsCoordinates(1, 3) = -r / norm;  mesh.Cell0DsCoordinates(2, 3) = -r / norm;
	mesh.Cell0DsCoordinates(0, 4) = -r / norm;  mesh.Cell0DsCoordinates(1, 4) = r / norm;  mesh.Cell0DsCoordinates(2, 4) = r / norm;
	mesh.Cell0DsCoordinates(0, 5) = -r / norm;  mesh.Cell0DsCoordinates(1, 5) = r / norm;  mesh.Cell0DsCoordinates(2, 5) = -r / norm;
	mesh.Cell0DsCoordinates(0, 6) = -r / norm;  mesh.Cell0DsCoordinates(1, 6) = -r / norm;  mesh.Cell0DsCoordinates(2, 6) = r / norm;
	mesh.Cell0DsCoordinates(0, 7) = -r / norm;  mesh.Cell0DsCoordinates(1, 7) = -r / norm;  mesh.Cell0DsCoordinates(2, 7) = -r / norm;
	
	double norm1= sqrt(3.0+sqrt(5.0));
    mesh.Cell0DsCoordinates(0, 8) = 0;  mesh.Cell0DsCoordinates(1, 8) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 8) = phi / norm1;
	mesh.Cell0DsCoordinates(0, 9) = 0;  mesh.Cell0DsCoordinates(1, 9) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 9) = -phi / norm1;
	mesh.Cell0DsCoordinates(0, 10) = 0;  mesh.Cell0DsCoordinates(1, 10) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 10) = phi / norm1;
	mesh.Cell0DsCoordinates(0, 11) = 0;  mesh.Cell0DsCoordinates(1, 11) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 11) = -phi / norm1;
	mesh.Cell0DsCoordinates(0, 12) = phi / norm1;  mesh.Cell0DsCoordinates(1, 12) = 0;  mesh.Cell0DsCoordinates(2, 12) = (r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 13) = -phi / norm1;  mesh.Cell0DsCoordinates(1, 13) = 0;  mesh.Cell0DsCoordinates(2, 13) = (r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 14) = phi / norm1;  mesh.Cell0DsCoordinates(1, 14) = 0;  mesh.Cell0DsCoordinates(2, 14) = (-r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 15) = -phi / norm1;  mesh.Cell0DsCoordinates(1, 15) = 0;  mesh.Cell0DsCoordinates(2, 15) = (-r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 16) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 16) = phi / norm1;  mesh.Cell0DsCoordinates(2, 16) = 0;
	mesh.Cell0DsCoordinates(0, 17) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 17) = -phi / norm1;  mesh.Cell0DsCoordinates(2, 17) = 0;
	mesh.Cell0DsCoordinates(0, 18) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 18) = phi / norm1;  mesh.Cell0DsCoordinates(2, 18) = 0;
	mesh.Cell0DsCoordinates(0, 19) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 19) = -phi / norm1;  mesh.Cell0DsCoordinates(2, 19) = 0;

    mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    
    // LATI
    mesh.Cell1DsId.reserve(30);
    Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    mesh.Cell1DsExtrema.resize(30, 2);
    mesh.Cell1DsExtrema <<
        0, 1,  // Lato 0
        1, 2,  // Lato 1
        2, 3,  // Lato 2
        3, 4,  // Lato 3
        4, 0,  // Lato 4
        1, 5,  // Lato 5
        6, 5,  // Lato 6
        6, 7,  // Lato 7
        7, 2,  // Lato 8
        6, 10,  // Lato 9
        8, 9,  // Lato 10
        9, 10,  // Lato 11
        10, 11,  // Lato 12
        5, 12,  // Lato 13
        11, 12,  // Lato 14
        12, 13,  // Lato 15
        15, 11,  // Lato 16
        14, 15,  // Lato 17
        13, 14,  // Lato 18
        0, 13,  // Lato 19
        15, 16,  // Lato 20
        14, 18,  // Lato 21
        18, 4,  // Lato 22
        18, 17,  // Lato 23
        17, 16,  // Lato 24
        17, 19,  // Lato 25
        16, 9,  // Lato 26
        19, 8,  // Lato 27
        19, 3,  // Lato 28
        8, 7;  // Lato 29

    // FACCE 
    mesh.Cell2DsId.reserve(12);
    mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    
    mesh.Cell2DsVertices.reserve(12, 5);
    mesh.Cell2DsEdges.reserve(12, 5);

     mesh.Cell2DsVertices = {
        {0, 1, 2, 3, 4},  // Faccia 0
        {1, 5, 6, 7, 2},  // Faccia 1
        {2, 7, 8, 19, 3},  // Faccia 2
        {10, 6, 7, 8, 9},  // Faccia 3
        {9, 8, 19, 17, 16},  // Faccia 4
        {3, 19, 17, 18, 4},  // Faccia 5
        {16, 17, 18, 14, 15},  // Faccia 6
        {15, 14, 13, 12, 11},   // Faccia 7
        {0, 4, 18, 14, 13},  // Faccia 8
        {13, 12, 5, 1, 0},  // Faccia 9
        {11, 12, 5, 6, 10},  // Faccia 10
        {9, 16, 15, 11, 10},  // Faccia 11
    };

    mesh.Cell2DsEdges = {
        {0, 1, 2, 3, 4},  // Faccia 0
        {5, 6, 7, 8, 1},  // Faccia 1
        {8, 29, 27, 28, 2},  // Faccia 2
        {9, 7, 29, 10, 11},  // Faccia 3
        {10, 27, 25, 24, 26},  // Faccia 4
        {28, 25, 23, 22, 3},  // Faccia 5
        {24, 23, 21, 17, 20},  // Faccia 6
        {17, 18, 15, 14, 16},   // Faccia 7
        {4, 22, 21, 18, 19},  // Faccia 8
        {15, 13, 5, 0, 19},  // Faccia 9
        {14, 13, 6, 9, 12},  // Faccia 10
        {26, 20, 16, 12, 11},  // Faccia 11
    };

    // POLIEDRO
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(20);
    mesh.Cell3DsEdges.reserve(30);
    mesh.Cell3DsFaces.reserve(12);
    
    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5,6 ,7 ,8, 9, 
	    10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
	    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    
    return true;
}

bool generateIcosahedron(PolyhedralMesh& mesh) {
    // VERTICI
    double r = 1.0;
    double phi = (1.0 + sqrt(5.0)) / 2.0;

    mesh.Cell0DsCoordinates.reserve(12, 3);
    mesh.Cell0DsId.reserve(12);
	
	double norm = sqrt(10.0+2.0*sqrt(5.0))/2.0;
    mesh.Cell0DsCoordinates(0, 0) = 0.0; mesh.Cell0DsCoordinates(1, 0) = r / norm; mesh.Cell0DsCoordinates(2, 0) = phi / norm;
    mesh.Cell0DsCoordinates(0, 1) = 0.0; mesh.Cell0DsCoordinates(1, 1) = -r / norm; mesh.Cell0DsCoordinates(2, 1) = phi/ norm;
    mesh.Cell0DsCoordinates(0, 2) = 0.0; mesh.Cell0DsCoordinates(1, 2) = r/ norm; mesh.Cell0DsCoordinates(2, 2) = -phi/ norm;
    mesh.Cell0DsCoordinates(0, 3) = 0.0; mesh.Cell0DsCoordinates(1, 3) = -r/ norm; mesh.Cell0DsCoordinates(2, 3) = -phi/ norm;
    mesh.Cell0DsCoordinates(0, 4) = r/ norm; mesh.Cell0DsCoordinates(1, 4) = phi/ norm; mesh.Cell0DsCoordinates(2, 4) = 0.0;
    mesh.Cell0DsCoordinates(0, 5) = -r/ norm; mesh.Cell0DsCoordinates(1, 5) = phi/ norm; mesh.Cell0DsCoordinates(2, 5) = 0.0;
    mesh.Cell0DsCoordinates(0, 6) = r/ norm; mesh.Cell0DsCoordinates(1, 6) = -phi/ norm; mesh.Cell0DsCoordinates(2, 6) = 0.0;
    mesh.Cell0DsCoordinates(0, 7) = -r/ norm; mesh.Cell0DsCoordinates(1, 7) = -phi/ norm; mesh.Cell0DsCoordinates(2, 7) = 0.0;
    mesh.Cell0DsCoordinates(0, 8) = phi/ norm; mesh.Cell0DsCoordinates(1, 8) = 0.0; mesh.Cell0DsCoordinates(2, 8) = r/ norm;
    mesh.Cell0DsCoordinates(0, 9) = -phi/ norm; mesh.Cell0DsCoordinates(1, 9) = 0.0; mesh.Cell0DsCoordinates(2, 9) = r/ norm;
    mesh.Cell0DsCoordinates(0, 10) = phi/ norm; mesh.Cell0DsCoordinates(1, 10) = 0.0; mesh.Cell0DsCoordinates(2, 10) = -r/ norm;
    mesh.Cell0DsCoordinates(0, 11) = -phi/ norm; mesh.Cell0DsCoordinates(1, 11) = 0.0; mesh.Cell0DsCoordinates(2, 11) = -r/ norm;
/ norm
    mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    // LATI
    mesh.Cell1DsId.reserve(30);
    mesh.Cell1DsExtrema.resize(30, 2);
    mesh.Cell1DsExtrema <<
    	0, 1,  // Lato 0
        1, 2,  // Lato 1
        2, 0,  // Lato 2
        3, 0,  // Lato 3
        4, 0,  // Lato 4
        0, 5,  // Lato 5
        1, 5,  // Lato 6
        3, 2,  // Lato 7
        5, 4,  // Lato 8
        4, 3,  // Lato 9
        4, 7,  // Lato 10
        5, 7,  // Lato 11
        5, 6,  // Lato 12
        1, 6,  // Lato 13
        1, 10,  // Lato 14
        4, 8,  // Lato 15
        2, 10,  // Lato 16
        3, 8,  // Lato 17
        2, 9,  // Lato 18
        6, 10,  // Lato 19
        10, 9,  // Lato 20
        6, 7,  // Lato 21
        7, 8,  // Lato 22
        8, 9,  // Lato 23
        7, 11,  // Lato 24
        8, 11,  // Lato 25
        10, 11,  // Lato 26
        6, 11,  // Lato 27
        9, 11,  // Lato 28
        3, 9;  // Lato 29

	// FACCE
	mesh.Cell2DsId.reserve(20);
	mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
	
	mesh.Cell2DsVertices.reserve(20, 3);
	mesh.Cell2DsVertices = {
		{0, 1, 2},  // Faccia 0
        {0, 2, 3},  // Faccia 1
        {0, 3, 4},  // Faccia 2
        {0, 4, 5},  // Faccia 3
        {0, 5, 1},  // Faccia 4
        {5, 7, 4},  // Faccia 5
        {4, 7, 8},  // Faccia 6
        {4, 8, 3},   // Faccia 7
        {3, 8, 9},  // Faccia 8
        {9, 3, 2},  // Faccia 9
        {2, 10, 9},  // Faccia 10
        {2, 1, 10},  // Faccia 11
        {1, 6, 10},  // Faccia 12
        {1, 5, 6},  // Faccia 13
        {5, 6, 7},  // Faccia 14
        {6, 7, 11},  // Faccia 15
        {7, 11, 8},  // Faccia 16
        {8, 11, 9},  // Faccia 17
        {9, 10, 11},  // Faccia 18
        {6, 10, 11},   // Faccia 19
    };
    mesh.Cell2DsEdges.reserve(20,3);
	mesh.Cell2DsEdges = {
		{0, 1, 2},  // Faccia 0
        {2, 7, 3},  // Faccia 1
        {3, 9, 4},  // Faccia 2
        {4, 8, 5},  // Faccia 3
        {5, 6, 0},  // Faccia 4
        {11, 10, 8},  // Faccia 5
        {10, 22, 15},  // Faccia 6
        {15, 17, 9},   // Faccia 7
        {17, 23, 29},  // Faccia 8
        {29, 7, 18},  // Faccia 9
        {16, 20, 18},  // Faccia 10
        {1, 14, 16},  // Faccia 11
        {13, 19, 14},  // Faccia 12
        {6, 12, 13},  // Faccia 13
        {12, 21, 11},  // Faccia 14
        {21, 24, 27},  // Faccia 15
        {24, 25, 22},  // Faccia 16
        {25, 28, 23},  // Faccia 17
        {20, 26, 28},  // Faccia 18
        {19, 26, 27},   // Faccia 19
    };
    
    // POLIEDRO
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(12);
    mesh.Cell3DsEdges.reserve(30);
    mesh.Cell3DsFaces.reserve(20);
    
    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5,6 ,7 ,8, 9, 10, 11};
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
	    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    
    return true;
}


bool subdivideTriangle(PolyhedralMesh& mesh, unsigned int b, unsigned int c) {
    if ((b == 0 && c == 0) || (b != 0 && c != 0)) return false;

    std::ofstream output0D("Cell0Ds.txt");
    std::ofstream output1D("Cell1Ds.txt");
    std::ofstream output2D("Cell2Ds.txt");
    std::ofstream output3D("Cell3Ds.txt");

    if (!output0D || !output1D || !output2D || !output3D) return false;

    std::vector<Eigen::Vector2d> vertices;
    std::vector<std::array<int, 3>> triangles;

    for (const auto& face : mesh.Cell2DsVertices) {
        Eigen::Vector2d A = mesh.Cell0DsCoordinates.row(face[0]);
        Eigen::Vector2d B = mesh.Cell0DsCoordinates.row(face[1]);
        Eigen::Vector2d C = mesh.Cell0DsCoordinates.row(face[2]);

        std::vector<std::vector<int>> vertexGrid;
		
		if (c == 0){
        for (int i = 0; i <= b; ++i) {
            std::vector<int> row;
            Eigen::Vector2d start = ((double)i / b) * A + ((double)(b - i) / b) * C;
            Eigen::Vector2d end = ((double)i / b) * B + ((double)(b - i) / b) * C;
            for (int j = 0; j <= i; ++j) {
                Eigen::Vector2d point;
                if (i == 0) {
                    point = C;
                } else {
                    point = ((double)j / i) * end + ((double)(i - j) / i) * start;
                }
                vertices.push_back(point);
                row.push_back(vertices.size() - 1);
            }
            vertexGrid.push_back(row);
        }

        // Ora costruiamo i triangoli
        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < i; ++j) {
                int v1 = vertexGrid[i][j];
                int v2 = vertexGrid[i + 1][j];
                int v3 = vertexGrid[i + 1][j + 1];
                triangles.push_back({v1, v2, v3});

                int v4 = vertexGrid[i][j];
                int v5 = vertexGrid[i + 1][j + 1];
                int v6 = vertexGrid[i][j + 1];
                triangles.push_back({v4, v5, v6});
            }
            // Triangolo in fondo a sinistra
            int v1 = vertexGrid[i][i];
            int v2 = vertexGrid[i + 1][i];
            int v3 = vertexGrid[i + 1][i + 1];
            triangles.push_back({v1, v2, v3});
        }
		}
		else {
		for (int i = 0; i <= c; ++i) {
			std::vector<int> row;
			Eigen::Vector2d start = ((double)i / c) * A + ((double)(c - i) / c) * B;
			Eigen::Vector2d end   = ((double)i / c) * C + ((double)(c - i) / c) * B;

			for (int j = 0; j <= i; ++j) {
				Eigen::Vector2d point;
				if (i == 0) {
					point = B;
				} else {
					point = ((double)j / i) * end + ((double)(i - j) / i) * start;
				}
				vertices.push_back(point);
				row.push_back(vertices.size() - 1);
			}
			vertexGrid.push_back(row);
		}

		// Triangoli
		for (int i = 0; i < c; ++i) {
			for (int j = 0; j < i; ++j) {
				int v1 = vertexGrid[i][j];
				int v2 = vertexGrid[i + 1][j];
				int v3 = vertexGrid[i + 1][j + 1];
				triangles.push_back({v1, v2, v3});

				int v4 = vertexGrid[i][j];
				int v5 = vertexGrid[i + 1][j + 1];
				int v6 = vertexGrid[i][j + 1];
				triangles.push_back({v4, v5, v6});
			}
			// Triangolo finale a sinistra
			int v1 = vertexGrid[i][i];
			int v2 = vertexGrid[i + 1][i];
			int v3 = vertexGrid[i + 1][i + 1];
			triangles.push_back({v1, v2, v3});
		}
		}
    }

    // Scrittura su file output (esempio Cell0Ds)
    for (size_t i = 0; i < vertices.size(); ++i) {
        output0D << i << " " << vertices[i].x() << " " << vertices[i].y() << std::endl;
    }

    for (size_t i = 0; i < triangles.size(); ++i) {
        output2D << i << " " << triangles[i][0] << " " << triangles[i][1] << " " << triangles[i][2] << std::endl;
    }

    output0D.close();
    output1D.close();
    output2D.close();
    output3D.close();
    return true;
}

}