#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "UCDUtilities.hpp"
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
namespace PolyhedralLibrary
{
	void ExportParaview(const PolyhedralMesh& meshFinal){

    Gedim::UCDUtilities utilities;

/*    bool hasValidVertexMarkers = !meshFinal.Cell0DsMarker.empty() &&
                                 (meshFinal.Cell0DsMarker.size() == meshFinal.Cell0DsCoordinates.rows());

    bool hasValidEdgeMarkers = !meshFinal.Cell1DsMarker.empty() &&
                              (meshFinal.Cell1DsMarker.size() == meshFinal.Cell1DsExtrema.rows()); */
	
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

        utilities.ExportPoints("./Cell0Ds.inp",
                              meshFinal.Cell0DsCoordinates,
                              { vertexProperty },
                              Eigen::VectorXi());

        utilities.ExportSegments("./Cell1Ds.inp",
                                meshFinal.Cell0DsCoordinates,
                                meshFinal.Cell1DsExtrema.transpose(),
                                { vertexProperty },
                                { edgeProperty },
                                Eigen::VectorXi());

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
	
		cout << "\n--- Fine struttura ---" << endl;
	}
}