#include "MeshExport.hpp"
#include <fstream>
#include <sstream>

using namespace std;
namespace PolyhedralLibrary {

    void WriteCell0DCSV(const PolyhedralMesh& mesh, const string& Cell0D) {
        ofstream file(Cell0D);
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file " << Cell0D << " per la scrittura." << endl;
        return;
    }
        file << "Id,x,y,z\n";
        for (size_t i = 0; i < mesh.Cell0DsId.size(); ++i) {
            file << mesh.Cell0DsId[i] << ","
                 << mesh.Cell0DsCoordinates(0, i) << ","
                 << mesh.Cell0DsCoordinates(1, i) << ","
                 << mesh.Cell0DsCoordinates(2, i) << "\n";
        }
        file.close();
    }

    void WriteCell1DCSV(const PolyhedralMesh& mesh, const string& Cell1D) {
        ofstream file(Cell1D);
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file " << Cell1D << " per la scrittura." << endl;
        return;
    }
        file << "Id,Start,End\n";
        for (size_t i = 0; i < mesh.Cell1DsId.size(); ++i) {
            file << mesh.Cell1DsId[i] << ","
                 << mesh.Cell1DsExtrema(i, 0) << ","
                 << mesh.Cell1DsExtrema(i, 1) << "\n";
        }
        file.close();
    }

    void WriteCell2DCSV(const PolyhedralMesh& mesh, const string& Cell2D){
        std::ofstream file(Cell2D);
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file " << Cell2D << " per la scrittura." << endl;
        return;
    }
        file << "Id,Vertices,Edges\n";
        for (size_t i = 0; i < mesh.Cell2DsId.size(); ++i) {
            file << mesh.Cell2DsId[i] << ",";
            for (size_t j = 0; j < mesh.Cell2DsVertices[i].size(); ++j) {
                file << mesh.Cell2DsVertices[i][j];
                if (j != mesh.Cell2DsVertices[i].size() - 1) file << "-"; // separo i valori con i trattini
            }
            file << ",";
            for (size_t j = 0; j < mesh.Cell2DsEdges[i].size(); ++j) {
                file << mesh.Cell2DsEdges[i][j];
                if (j != mesh.Cell2DsEdges[i].size() - 1) file << "-";
            }
            file << "\n";
        }
        file.close();
    }

    void WriteCell3DCSV(const PolyhedralMesh& mesh, const string& Cell3D) {
        ofstream file(Cell3D);
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file " << Cell3D << " per la scrittura." << endl;
        return;
    }
        file << "Id,Vertices,Edges,Faces\n";
        for (size_t i = 0; i < mesh.Cell3DsId.size(); ++i) {
            file << mesh.Cell3DsId[i] << ",";
            file << mesh.Cell3DsVertices[i] << ",";
            file << mesh.Cell3DsEdges[i] << ",";
            file << mesh.Cell3DsFaces[i] << "\n";
        }
        file.close();
    }

}
