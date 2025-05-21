#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>

using namespace std;
namespace PolyhedralLibrary {

    void WriteCell0Ds(PolyhedralMesh& mesh) {
        ofstream file("Cell0Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell0Ds.txt per la scrittura." << endl;
        return;
    	}
    	unsigned int maxFlag = numeric_limits<unsigned int>::max();
        file << "Id;x;y;z;duplicated\n";
        for (size_t i = 0; i < mesh.Cell0DsId.size(); ++i) {
            file << mesh.Cell0DsId[i] << ";"
                 << mesh.Cell0DsCoordinates(0, i) << ";"
                 << mesh.Cell0DsCoordinates(1, i) << ";"
                 << mesh.Cell0DsCoordinates(2, i) << ";";
            if ( mesh.Cell0DsFlag[i][0] == maxFlag) {
	            file << "false \n";
	        }
	        else {
		        file << "true \n";
        	}
    	}
    file.close();
    }

    void WriteCell1Ds(PolyhedralMesh& mesh) {
        ofstream file("Cell1Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell1Ds.txt per la scrittura." << endl;
        return;
    	}
    	unsigned int maxFlag = numeric_limits<unsigned int>::max();
        file << "Id;Start;End;duplicated\n";
        for (size_t i = 0; i < mesh.Cell1DsId.size(); ++i) {
            file << mesh.Cell1DsId[i] << ";"
                 << mesh.Cell1DsExtrema(i, 0) << ";"
                 << mesh.Cell1DsExtrema(i, 1) << ";";
            if (mesh.Cell1DsFlag[i] == maxFlag){
	            file << "false \n";
	        }
	        else {
		        file << "true \n";
		    }
        }
        file.close();
    }

    void WriteCell2Ds(PolyhedralMesh& mesh){
        std::ofstream file("Cell2Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell2Ds.txt per la scrittura." << endl;
        return;
    }
        file << "Id;NumVertices;Vertices;NumEdges;Edges\n";
        for (size_t i = 0; i < mesh.Cell2DsId.size(); ++i) {
            file << mesh.Cell2DsId[i] << ";" << "3;";
            for (size_t j = 0; j < mesh.Cell2DsVertices[i].size(); ++j) {
                file << mesh.Cell2DsVertices[i][j] << ";";
            }
            file << "3;";
            for (size_t j = 0; j < mesh.Cell2DsEdges[i].size(); ++j) {
                file << mesh.Cell2DsEdges[i][j] << ";";
            }
            file << "\n";
        }
        file.close();
    }

    void WriteCell3Ds(PolyhedralMesh& mesh) {
        ofstream file("Cell3Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell3Ds.txt per la scrittura." << endl;
        return;
    }
        file << "Id;Vertices;Edges;Faces\n";
        for (size_t i = 0; i < mesh.Cell3DsId.size(); ++i) {
            file << mesh.Cell3DsId[i] << ";";
            file << mesh.Cell3DsVertices[i] << ";";
            file << mesh.Cell3DsEdges[i] << ";";
            file << mesh.Cell3DsFaces[i] << "\n";
        }
        file.close();
    }

}
