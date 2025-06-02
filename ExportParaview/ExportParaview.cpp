#include "UCDUtilities.hpp"
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <Eigen/Dense>

using namespace std;
namespace PolyhedralLibrary
{
	void ExportParaview(const PolyhedralMesh& meshFinal){

    Gedim::UCDUtilities utilities;
	
	Eigen::VectorXi vertexIds(meshFinal.Cell0DsId.size());
	for (size_t i = 0; i < meshFinal.Cell0DsId.size(); ++i)
		vertexIds[i] = static_cast<int>(meshFinal.Cell0DsId[i]);

	Eigen::VectorXi edgeIds(meshFinal.Cell1DsId.size());
	for (size_t i = 0; i < meshFinal.Cell1DsId.size(); ++i)
		edgeIds[i] = static_cast<int>(meshFinal.Cell1DsId[i]);
	
	bool hasValidVertexMarkers = !meshFinal.Cell0DsMarker.empty();

    bool hasValidEdgeMarkers = !meshFinal.Cell1DsMarker.empty();

    if (hasValidVertexMarkers && hasValidEdgeMarkers) {
        // Creo e assegno le proprietà scalari
        vector<double> vertexData;
        vertexData.reserve(meshFinal.Cell0DsMarker.size());
        for (auto m : meshFinal.Cell0DsMarker)
            vertexData.push_back(static_cast<double>(m));

        Gedim::UCDProperty<double> vertexProperty;
        vertexProperty.Label = "PathVertex";
        vertexProperty.UnitLabel = "-";
        vertexProperty.NumComponents = 1;
		
        vertexProperty.Data = vertexData.data();
        vertexProperty.Size = vertexData.size();

        vector<double> edgeData;
        edgeData.reserve(meshFinal.Cell1DsMarker.size());
        for (auto m : meshFinal.Cell1DsMarker)
            edgeData.push_back(static_cast<double>(m));

        Gedim::UCDProperty<double> edgeProperty;
        edgeProperty.Label = "PathEdge";
        edgeProperty.UnitLabel = "-";
        edgeProperty.NumComponents = 1;
        edgeProperty.Data = edgeData.data();
        edgeProperty.Size = edgeData.size();
		
		/*utilities.ExportPoints("./Cell0Ds.inp",
                              meshFinal.Cell0DsCoordinates,
                              { vertexProperty },
                              Eigen::VectorXi());*/

        utilities.ExportPoints("./Cell0Ds.inp",
                              meshFinal.Cell0DsCoordinates,
                              { vertexProperty },
                              vertexIds);

        utilities.ExportSegments("./Cell1Ds.inp",
                                meshFinal.Cell0DsCoordinates,
                                meshFinal.Cell1DsExtrema.transpose(),
                                { vertexProperty },
                                { edgeProperty },
                                edgeIds);

    } else {
        // Export senza proprietà scalari
        utilities.ExportPoints("./Cell0Ds.inp",
                              meshFinal.Cell0DsCoordinates,
                              {},  // nessuna proprietà
                              Eigen::VectorXi());

        utilities.ExportSegments("./Cell1Ds.inp",
                                meshFinal.Cell0DsCoordinates,
                                meshFinal.Cell1DsExtrema.transpose(),
                                {}, // nessuna proprietà vertici
                                {}, // nessuna proprietà segmenti
                                Eigen::VectorXi());
    }
	
	}
	
	void printMeshTriangulated(const PolyhedralMesh& mesh) {
		
		// VERTICI
		cout << "Cell0DsId: "; 
		for (auto id : mesh.Cell0DsId) cout << id << " ";
		cout << "\nCell0DsCoordinates (per colonne):" << endl;
		for (int j = 0; j < mesh.Cell0DsCoordinates.cols(); ++j) {
			cout << "Colonna " << j << ": ";
			for (int i = 0; i < mesh.Cell0DsCoordinates.rows(); ++i) {
				cout << mesh.Cell0DsCoordinates(i, j) << " ";
			}
			cout << endl;
		}
		cout << "Cell0DsFlag:" << endl;
		for (const auto& row : mesh.Cell0DsFlag) {
			for (auto v : row) cout << v << " ";
			cout << endl;
		}
		
		// LATI
		cout << "Cell1DsId: "; 
		for (auto id : mesh.Cell1DsId) cout << id << " ";
		cout << "\nCell1DsExtrema (per righe):" << endl;
		for (int i = 0; i < mesh.Cell1DsExtrema.rows(); ++i) {
			for (int j = 0; j < mesh.Cell1DsExtrema.cols(); ++j) {
				cout << mesh.Cell1DsExtrema(i, j) << " ";
			}
			cout << endl;
		}
		cout << "Cell1DsFlag:" << endl;
		for (const auto& row : mesh.Cell1DsFlag) {
			cout << row << " ";
			cout << endl;
		}
		
		// FACCE
		cout << "Cell2DsId: "; 
		for (auto id : mesh.Cell2DsId) cout << id << " ";
		
		cout << "\nCell2DsVertices:" << endl;
		for (const auto& row : mesh.Cell2DsVertices) {
			for (auto v : row) cout << v << " ";
			cout << endl;
		}
		cout << "Cell2DsEdges:" << endl;
		for (const auto& row : mesh.Cell2DsEdges) {
			for (auto v : row) cout << v << " ";
			cout << endl;
		}
		
		// POLIEDRO
		cout << "Cell3DsId: " << mesh.Cell3DsId << endl; 
		cout << "Numero di vertici: " << mesh.NumCells0Ds << endl;
		cout << "Numero di lati: " << mesh.NumCells1Ds << endl;
		cout << "Numero di facce: " << mesh.NumCells2Ds << endl;
		
		cout << "Cell3DsVertices: "; 
		for (auto id : mesh.Cell3DsVertices) cout << id << " ";
		cout << "Cell3DsEdges: "; 
		for (auto id : mesh.Cell3DsEdges) cout << id << " ";
		cout << "Cell3DsFaces: "; 
		for (auto id : mesh.Cell3DsFaces) cout << id << " ";
	
		cout << "\n--- Fine struttura ---" << endl;
	}
}