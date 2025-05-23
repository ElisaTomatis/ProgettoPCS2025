#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

namespace PolyhedralLibrary
{
	vector<int> ComputePolyhedronVEF(int q, int b, int c);
	
	vector<int> CalculateDuplicated(int q, int b, int c, const vector<int>& dimension);
	
	void RemoveDuplicatedVertices(PolyhedralLibrary::PolyhedralMesh& meshTriangulated);
	
	void RemoveDuplicatedEdges(PolyhedralLibrary::PolyhedralMesh& meshTriangulated);
	
	void generateTetrahedron(PolyhedralLibrary::PolyhedralMesh& mesh);
	
	void generateCube(PolyhedralLibrary::PolyhedralMesh& mesh);
	
	void generateOctahedron(PolyhedralLibrary::PolyhedralMesh& mesh);
	
	void generateDodecahedron(PolyhedralLibrary::PolyhedralMesh& mesh);
	
	void generateIcosahedron(PolyhedralLibrary::PolyhedralMesh& mesh);
	
	/**
	 * Triangola una faccia del poliedro e salva il risultato nella mesh triangolata.
	 * 
	 * mesh Mesh di input contenente le facce da triangolare.
	 * meshTriangulated Mesh risultante con le facce triangolate.
	 * ID della faccia da triangolare.
	 * cellID ID della cella a cui appartiene la faccia.
	 * dimension Informazioni di dimensione, tipicamente [numVertices, numEdges, numFaces].
	 */
	 
	 
	void triangulateAndStore(PolyhedralLibrary::PolyhedralMesh& mesh, PolyhedralLibrary::PolyhedralMesh& meshTriangulated,
							  unsigned int b, unsigned int c, const vector<int>& dimension);
	
	
	/**
	 * Aggiunge un lato alla mesh triangolata se non già presente.
	 * 
	 * a Primo estremo del lato.
	 * b Secondo estremo del lato.
	 * meshTriangulated Mesh in cui aggiungere il lato.
	 * edgeID ID globale del nuovo lato (incrementato se un lato viene aggiunto).
	 * triangleID ID del triangolo a cui il lato appartiene.
	 */
	void FindAddEdge(unsigned int a, unsigned int b, PolyhedralLibrary::PolyhedralMesh& meshTriangulated,
					 unsigned int& edgeID, unsigned int triangleID);
					 
	
	/**
	* Esporta la mesh triangolata a paraview.
	*/
	void ExportParaview(const PolyhedralLibrary::PolyhedralMesh& meshFinal);
	
	void printMeshTriangulated(const PolyhedralLibrary::PolyhedralMesh& meshTriangulated);
	
	void WriteCell0Ds(PolyhedralLibrary::PolyhedralMesh& mesh);
    void WriteCell1Ds(PolyhedralLibrary::PolyhedralMesh& mesh);
    void WriteCell2Ds(PolyhedralLibrary::PolyhedralMesh& mesh);
    void WriteCell3Ds(PolyhedralLibrary::PolyhedralMesh& mesh);
    
    void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual);
    map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTrinagulated);
	
}