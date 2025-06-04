#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>


using namespace PolyhedralLibrary;
using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	// DIMENSION
	
	// Inverte i valori di p e q
	void invertiValori(int& p, int& q);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato con la triangolazione I
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	vector<int> ComputePolyhedronVEF(int q, int b, int c);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato (considerati i duplicati)
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato (senza duplicati)
	vector<int> CalculateDuplicated(int q, int b, int c, const vector<int>& dimension);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato con la triangolazione II
	// q, b : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	vector<int> CalculateDimension2(int b, int q);
	
	// Assegna un flag ai vertici che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct PolyhedralMesh triangolata
	void RemoveDuplicatedVertices(PolyhedralMesh& meshTriangulated);
	
	// Assegna un flag ai lati che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct PolyhedralMesh triangolata
	void RemoveDuplicatedEdges(PolyhedralMesh& meshTriangulated);
	
	// Crea una nuova mesh senza duplicati
	// meshTriangulated : una struct PolyhedralMesh (con i duplicati)
	// meshFinal : una struct PolyhedralMesh (senza i duplicati)
	// dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato (senza duplicati)
	void NewMesh(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension);
	
	// ----------------------------------------------------------------------------------------------------------- // 
	
	// POLYHEDRA 
	
	// Riempiono la struct PolyhedralMesh con i dati dei poliedri non triangolati
	// mesh: una struct PolyhedralMesh
	void generateTetrahedron(PolyhedralMesh& mesh);
	void generateCube(PolyhedralMesh& mesh);
	void generateOctahedron(PolyhedralMesh& mesh);
	void generateDodecahedron(PolyhedralMesh& mesh);
	void generateIcosahedron(PolyhedralMesh& mesh);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// TRIANGULATION
	
	// Racchiude le funzioni per la triangolazione di classe I se p=3
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// mesh: una struct PolyhedralMesh
	PolyhedralMesh Triangulation(int q, int b, int c, PolyhedralMesh& mesh);
	
	// Racchiude le funzioni per la triangolazione di classe I se q=3
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// mesh: una struct PolyhedralMesh
	PolyhedralMesh TriangulationDual(int q, int b, int c, PolyhedralMesh& mesh);
	
	// Racchiude le funzioni per la triangolazione di classe II se p=3
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// mesh: una struct PolyhedralMesh
	PolyhedralMesh Triangulation2(int q, int b, int c, PolyhedralMesh& mesh);
	
	// Racchiude le funzioni per la triangolazione di classe II se q=3
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// mesh: una struct PolyhedralMesh
	PolyhedralMesh Triangulation2Dual(int q, int b, int c, PolyhedralMesh& mesh);
	
	// Riempie le Celle3d dopo la triangolazione
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato
	void PopulateCell3D(PolyhedralMesh& meshTriangulated, const vector<int>& dimension);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// TRIANGULATION1
	
	// Triangola il poliedro di classe I
	// mesh : una struct PolyhedralMesh, quella di partenza non triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// b,c  : parametri passati dall'utente che identificano il poliedro
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
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// TRIANGULATION2
					 
	// Triangola il poliedro di classe II
	 // mesh : una struct PolyhedralMesh, quella di partenza non triangolata
	 // meshTriangulated : una struct PolyhedralMesh, quella triangolata
	 // b,c : parametri passati dall'utente che identificano il poliedro
	 // dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato
	void triangulateAndStore2(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated, const vector<int>& dimension);
							  
	// Aggiunge un vertice alla mesh triangolata se non è già presente e restituisce il suo id
	// coord : coordinate del punto che vorremmo aggiungere
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// k1 : prossimo id per i vertici disponibile 
	unsigned int FindAddVertice(const Vector3d& coord, PolyhedralMesh& meshTriangulated, unsigned int& k1);
	
	// Aggiunge un lato alla mesh triangolata se non è già presente
	// a : id del primo estremo (vertice) del lato
	// b : id del secondo estremo (vertice) del lato
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// k2 : prossimo id per i lati disponibile
	unsigned int FindAddEdge2(unsigned int a, unsigned int b, PolyhedralMesh& meshTriangulated, unsigned int& k2);
	
	// A partire da un lato e una faccia che stiamo considerando, trova l'altra faccia adiacente
	// al lato e restituisce il suo baricentro
	// meshTriangulated : una struct PolyhedralMesh
	// edgeId : id del lato
	// currentFacdeId : faccia che stiamo considerando
	Vector3d FindNearBarycenter(const PolyhedralMesh& meshTriangulated, unsigned int edgeId, unsigned int currentFaceId);
	
	// Aggiunge una faccia alla mesh triangolata se non è già presente
	// new_face_vertices : id dei vertici della faccia che vorremmo aggiungere
	// new_face_edges : id dei lati della faccia che vorremmo aggiungere
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// k3 : prossimo id per le facce disponibile
	void FindAddFace(const vector<unsigned int>& new_face_vertices, const vector<unsigned int>& new_face_edges, PolyhedralMesh& meshTriangulated, unsigned int& k3);
	
	// Ruota la sequenza degli id dei lati in modo che inizi con il più piccolo (8,5,10) --> (5,10,8)
	// current_edges : id dei lati che vogliamo riordinare
	vector<unsigned int> get_cyclic_normalized(const vector<unsigned int>& current_edges);
	
	// Normalizza l'ordine combinando la normalizzazione ciclica con un confronto delle forme invertite (0,1,2) e (0,2,1) sono uguali
	// face_edges : id dei lati che vogliamo riordinare
	vector<unsigned int> NormalizeFaceEdges(const vector<unsigned int>& face_edges);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// DUAL
	
	// Calcola il duale di un poliedro
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// meshDual : una struct PolyhedralMesh, quella duale
	void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual);
	
	// Calcola il baricentro di una faccia
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// faceId : id della faccia di cui vogliamo calcolare il baricentro
	Vector3d getFaceBarycenter(const PolyhedralMesh& meshTriangulated, unsigned int faceId);
	
	// Crea una mappa che associa ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
    map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTriangulated);
	
	// Crea una mappa che associa ad ogni vertice originale l'elenco di tutte le facce che contengono quel vertice
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	map<unsigned int, vector<unsigned int>> buildVertexToFacesMap(const PolyhedralMesh& meshTriangulated);
	
	// Crea una mappa che associa ad ogni vertice originale l'elenco di tutti gli spigoli che vi incidono
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	map<unsigned int, vector<unsigned int>> buildVertexToEdgesMap(const PolyhedralMesh& meshTriangulated);
	
	// Proietta i vertici del poliedro sulla sfera unitaria
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void ProjectMeshToUnitSphere(PolyhedralMesh& meshTriangulated);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// CAMMINO MINIMO
	 
	struct ShortestPathResult {
		unsigned int numEdges;
		double totalLength;
		vector<bool> verticesInPath;
		vector<bool> edgesInPath;

		// Costruttore che inizializza i vettori con la dimensione corretta
		ShortestPathResult(unsigned int nEdges = 0, double len = 0.0,
                       unsigned int numV = 0, unsigned int numE = 0)
        : numEdges(nEdges), totalLength(len),
          verticesInPath(numV, false), edgesInPath(numE, false)
		{}
	};
	
	double calculateDistanceById(const PolyhedralMesh& mesh, const map<unsigned int, unsigned int>& vertexIdToIndexMap, unsigned int id1, unsigned int id2);
	MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh);
	
	ShortestPathResult findShortestPathDijkstra(
	    PolyhedralMesh& mesh,
	    const MatrixXi& adjMatrix,
	    unsigned int startVertexId_real,
	    unsigned int endVertexId_real
	);

	// pair<unsigned int, double> findShortestPathBFS(PolyhedralLibrary::PolyhedralMesh& mesh, const MatrixXi& adjMatrix, unsigned int startVertexId_real, unsigned int endVertexId_real, vector<bool>& isVertexInShortestPath, vector<bool>& isEdgeInShortestPath);
	
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// EXPORT PARAVIEW
	
	// Esporta la mesh triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void ExportParaview(const PolyhedralMesh& meshTriangulated);
	
	// Stampa la mesh triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void printMeshTriangulated(const PolyhedralMesh& meshTriangulated);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// MESH EXPORT

	// Scrivono sui file TXT
	// mesh: una struct PolyhedralMesh
	void WriteCell0Ds(const PolyhedralMesh& mesh);
    void WriteCell1Ds(const PolyhedralMesh& mesh);
    void WriteCell2Ds(const PolyhedralMesh& mesh);
    void WriteCell3Ds(const PolyhedralMesh& mesh);
    
    // ----------------------------------------------------------------------------------------------------------- //

}