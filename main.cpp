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
    int q = 5; //atoi(argv[2]);
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
    } else if (p == 3 && q == 4) {
        PolyhedralLibrary::generateOctahedron(mesh);
    } else if (p == 3 && q == 5) {
        PolyhedralLibrary::generateIcosahedron(mesh);
    } else {
        cerr << "Errore: combinazione p=" << p << " e q=" << q << " non supportata.\n";
        return 1;
    }

    vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
    vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
    PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
    PolyhedralLibrary::ExportParaview(meshTriangulated);

    return 0;
}

	
	