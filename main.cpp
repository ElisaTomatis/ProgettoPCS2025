#include <iostream>
#include <cstdlib> 
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main(int argc, char *argv[]) {
    if (argc != 5) { // argc conterrà il numero totale di argomenti, incluso il nome del programma stesso
        std::cerr << "Uso: " << argv[0] << " p q b c\n";
        return 1;
    }
    
    // argv[0] = nome del programma
    int p = stoi(argv[1]);
    int q = stoi(argv[2]);
    int b = stoi(argv[3]);
    int c = stoi(argv[4]);

    cout << "Hai inserito: " << p << " " << q << " " << b << " " << c << "\n";

    if ((p != 3 && p != 4 && p != 5) || (q != 3 && q != 4 && q != 5)) {
        cerr << "Errore: p e q devono essere 3, 4 o 5.\n";
        return 1;
    }

    PolyhedralLibrary::PolyhedralMesh mesh;
    PolyhedralLibrary::PolyhedralMesh meshTriangulated;
    PolyhedralLibrary::PolyhedralMesh meshFinal;
    PolyhedralLibrary::PolyhedralMesh meshDual;
    
    if (p == 3 && q == 3) {
        PolyhedralLibrary::generateTetrahedron(mesh);
        vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
		PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
    	PolyhedralLibrary::PopulateCell3D(meshTriangulated, dimension);
    	//printMeshTriangulated(meshTriangulated);
    	PolyhedralLibrary::ExportParaview(meshTriangulated);
		
    } else if (p == 3 && q != 3){
	    if (q == 4){
		    PolyhedralLibrary::generateOctahedron(mesh);
		} else {
			PolyhedralLibrary::generateIcosahedron(mesh);
		}
		vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
		PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
    	PolyhedralLibrary::PopulateCell3D(meshTriangulated, dimension);
    	//printMeshTriangulated(meshTriangulated);
    	PolyhedralLibrary::ExportParaview(meshTriangulated);

	} else if (q == 3 && p!= 3) {
		if (p == 4){
			PolyhedralLibrary::generateOctahedron(mesh);
		} else {
			PolyhedralLibrary::generateIcosahedron(mesh);
		}
		PolyhedralLibrary::invertiValori(p, q);
		vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
		PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
    	PolyhedralLibrary::CalculateDual(meshTriangulated, meshDual);
    	PolyhedralLibrary::PopulateCell3D(meshTriangulated, dimension);
    	//printMeshTriangulated(meshTriangulated);
    	PolyhedralLibrary::ExportParaview(meshDual);
		
	} else {
		cerr << "Errore: combinazione p =" << p << " e q =" << q << " non supportata.\n";
        return 1;
    }

    // Scrittura su TXT
	PolyhedralLibrary::WriteCell0Ds(meshTriangulated);
	PolyhedralLibrary::WriteCell1Ds(meshTriangulated);
	PolyhedralLibrary::WriteCell2Ds(meshTriangulated);
	PolyhedralLibrary::WriteCell3Ds(meshTriangulated);

    return 0;
}

	
	