#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "UCDUtilities.hpp"
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

namespace PolyhedralLibrary
{
	void ExportParaview(const PolyhedralMesh& meshFinal){
		
		Gedim::UCDUtilities utilities;
		{
			utilities.ExportPoints("./Cell0Ds.inp",
									 meshFinal.Cell0DsCoordinates,
									 {},
									 {});
	
			utilities.ExportSegments("./Cell1Ds.inp",
									 meshFinal.Cell0DsCoordinates,
									 meshFinal.Cell1DsExtrema.transpose(),
									 {},
									 {});
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