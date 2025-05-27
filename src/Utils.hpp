#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

using namespace PolyhedralLibrary;
namespace PolyhedralLibrary
{
	// Inverte i valori di p e q
	void invertiValori(int& p, int& q);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	vector<int> ComputePolyhedronVEF(int q, int b, int c);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato (considerati i duplicati)
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato
	vector<int> CalculateDuplicated(int q, int b, int c, const vector<int>& dimension);
	
	// Assegna un flag ai vertici che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct PolyhedralMesh
	void RemoveDuplicatedVertices(PolyhedralMesh& meshTriangulated);
	
	// Assegna un flag ai lati che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct PolyhedralMesh
	void RemoveDuplicatedEdges(PolyhedralMesh& meshTriangulated);
	
	// Riempiono la struct PolyhedralMesh con i dati dei poliedri non traingolati
	// mesh: una struct PolyhedralMesh
	void generateTetrahedron(PolyhedralMesh& mesh);
	void generateCube(PolyhedralMesh& mesh);
	void generateOctahedron(PolyhedralMesh& mesh);
	void generateDodecahedron(PolyhedralMesh& mesh);
	void generateIcosahedron(PolyhedralMesh& mesh);
	
	 // Triangola il poliedro
	 // mesh : una struct PolyhedralMesh, quella di partenza non triangolata
	 // meshTriangulated : una struct PolyhedralMesh, quella triangolata
	 // b,c : parametri passati dall'utente che identificano il poliedro
	 // dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato
	void triangulateAndStore(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated,
							  unsigned int b, unsigned int c, const vector<int>& dimension);
	
	 // Aggiunge un lato alla mesh triangolata se non è già presente
	 // a : id del primo estremo (vertice) del lato
	 // b : id del secondo estremo (vertice) del lato
	 // meshTriangulated : una struct PolyhedralMesh, quella triangolata
	 // edgeId : id del nuovo lato
	 // triangleId : id della faccia a cui appartiene
	void FindAddEdge(unsigned int a, unsigned int b, PolyhedralMesh& meshTriangulated,
					 unsigned int& edgeID, unsigned int triangleID);
	
	// Calcola il duale di un poliedro
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// meshDual : una struct PolyhedralMesh, quella duale
	void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual);
	
	// Crea una mappa che associa ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
    map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTriangulated);
	
	// Riempie le Celle3d dopo la triangolazione
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato
	void PopulateCell3D(PolyhedralMesh& meshTriangulated, const vector<int>& dimension);
	
	// Esporta la mesh triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void ExportParaview(const PolyhedralMesh& meshTriangulated);
	
	// Stampa la mesh triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void printMeshTriangulated(const PolyhedralMesh& meshTriangulated);
	
	// Proietta i vertici del poliedro sulla sfera unitaria
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void ProjectMeshToUnitSphere(PolyhedralMesh& mesh);
	
	// Scrivono sui file TXT
	// mesh: una struct PolyhedralMesh
	void WriteCell0Ds(const PolyhedralMesh& mesh);
    void WriteCell1Ds(const PolyhedralMesh& mesh);
    void WriteCell2Ds(const PolyhedralMesh& mesh);
    void WriteCell3Ds(const PolyhedralMesh& mesh);

	//double distance = calculateDistanceById(const PolyhedralMesh& mesh, const map<unsigned int, unsigned int>& vertexIdToIndexMap, unsigned int id1, unsigned int id2);
	//pair<unsigned int, double> path = findShortestPathBFS(const PolyhedralMesh& mesh, const MatrixXi& adjMatrix, unsigned int startVertexId_real, unsigned int endVertexId_real, vector<bool>& isVertexInShortestPath, vector<bool>& isEdgeInShortestPath);
	
}