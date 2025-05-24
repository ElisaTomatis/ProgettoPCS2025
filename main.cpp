#include <iostream>
#include <cstdlib> 
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main(int argc, char *argv[]) {
    //if (argc != 5) {
    //    std::cerr << "Uso: " << argv[0] << " p q b c\n";
    //    return 1;
    //}

    int p = 3; //atoi(argv[1]);
    int q = 3; //atoi(argv[2]);
    int b = 2; //atoi(argv[3]);
    int c = 0; //atoi(argv[4]);

    cout << "Hai inserito: " << p << " " << q << " " << b << " " << c << "\n";

    if ((p != 3 && p != 4 && p != 5) || (q != 3 && q != 4 && q != 5)) {
        cerr << "Errore: p e q devono essere 3, 4 o 5.\n";
        return 1;
    }

    PolyhedralLibrary::PolyhedralMesh mesh;
    PolyhedralLibrary::PolyhedralMesh meshTriangulated;
    
    if (p == 3 && q == 3) {
        PolyhedralLibrary::generateTetrahedron(mesh);
        vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
		PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
		
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

	} else if (q == 3 && p!= 3) {
		if (p == 4){
			PolyhedralLibrary::generateOctahedron(mesh);
		} else {
			PolyhedralLibrary::generateIcosahedron(mesh);
		}
		vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
		PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
		
	} else {
		cerr << "Errore: combinazione p =" << p << " e q =" << q << " non supportata.\n";
        return 1;
    }
    
    PolyhedralLibrary::ExportParaview(meshTriangulated);

    //PolyhedralLibrary::ExportParaview(mesh);
    PolyhedralLibrary::printMeshTriangulated(meshTriangulated);
    
    //PolyhedralLibrary::ExportParaview(meshTriangulated);
    
    // Scrittura su TXT
	//PolyhedralLibrary::WriteCell0Ds(meshTriangulated);
	//PolyhedralLibrary::WriteCell1Ds(meshTriangulated);
	//PolyhedralLibrary::WriteCell2Ds(meshTriangulated);
	//PolyhedralLibrary::WriteCell3Ds(meshTriangulated);

    return 0;
}

	
	