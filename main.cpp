#include <iostream>
#include <cstdlib> 
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main(int argc, char *argv[]) {
	/*
	
    // Definizione delle variabili per il cammino minimo, inizializzate a valori non validi
    unsigned int startVertexId = 0;
    unsigned int endVertexId = 0;
    bool calculatePath = false; // Flag per indicare se calcolare il cammino minimo

    if (argc == 5) { // Caso: solo p, q, b, c
        cout << "Modalità: Generazione mesh.\n";
        // Nessun cammino minimo da calcolare
        calculatePath = false;

    } else if (argc == 7) { // Caso: p, q, b, c, startVertexId, endVertexId
        cout << "Modalità: Generazione mesh e calcolo cammino minimo.\n";
        startVertexId = stoul(argv[5]);
        endVertexId = stoul(argv[6]);
        calculatePath = true;

    } else { // Errore: numero di argomenti non valido
        std::cerr << "Uso:\n";
        std::cerr << "  " << argv[0] << " p q b c\n";
        std::cerr << "  " << argv[0] << " p q b c startVertexId endVertexId\n";
        return 1;
    }
    
    // Parsing dei primi 4 argomenti, sempre presenti
    int p = stoi(argv[1]);
    int q = stoi(argv[2]);
    int b = stoi(argv[3]);
    int c = stoi(argv[4]);

    cout << "Hai inserito: p=" << p << ", q=" << q << ", b=" << b << ", c=" << c << "\n";
    if (calculatePath) {
        cout << "Cammino minimo da vertice ID: " << startVertexId << " a vertice ID: " << endVertexId << "\n";
    }

    if ((p != 3 && p != 4 && p != 5) || (q != 3 && q != 4 && q != 5)) {
        cerr << "Errore: p e q devono essere 3, 4 o 5.\n";
        return 1;
    }

    PolyhedralLibrary::PolyhedralMesh mesh;
    PolyhedralLibrary::PolyhedralMesh meshTriangulated;
    PolyhedralLibrary::PolyhedralMesh meshFinal;
    PolyhedralLibrary::PolyhedralMesh meshDual;
    
    if (p == 3 && q == 3) {
        PolyhedralLibrary::generateTetrahedron(mesh);
        vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
		PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
    	PolyhedralLibrary::PopulateCell3D(meshTriangulated, dimension);
    	//printMeshTriangulated(meshTriangulated);
    	// PolyhedralLibrary::ExportParaview(meshTriangulated);
		
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
    	PolyhedralLibrary::PopulateCell3D(meshTriangulated, dimension);
    	//printMeshTriangulated(meshTriangulated);
    	// PolyhedralLibrary::ExportParaview(meshTriangulated);

	} else if (q == 3 && p!= 3) {
		if (p == 4){
			PolyhedralLibrary::generateOctahedron(mesh);
		} else {
			PolyhedralLibrary::generateIcosahedron(mesh);
		}
		PolyhedralLibrary::invertiValori(p, q);
		vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
		PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
    	PolyhedralLibrary::CalculateDual(meshTriangulated, meshDual);
    	PolyhedralLibrary::PopulateCell3D(meshTriangulated, dimension);
    	//printMeshTriangulated(meshTriangulated);
    	// PolyhedralLibrary::ExportParaview(meshDual);
		
	} else {
		cerr << "Errore: combinazione p =" << p << " e q =" << q << " non supportata.\n";
        return 1;
    }
	
	// --- Selettore della mesh target per l'esportazione e il cammino minimo ---
    PolyhedralMesh* targetMeshPtr = nullptr;
    if (q == 3 && p != 3) {
        targetMeshPtr = &meshDual;
        if (calculatePath) cout << "Calcolo cammino minimo sulla mesh Duale.\n";
    } else {
        targetMeshPtr = &meshTriangulated;
        if (calculatePath) cout << "Calcolo cammino minimo sulla mesh Triangolata.\n";
    }

    if (targetMeshPtr->Cell0DsId.empty()) {
        cerr << "Errore: La mesh target (triangolata o duale) è vuota. Impossibile proseguire.\n";
        return 1;
    }
	
	// --- Calcolo del cammino minimo BFS (condizionale) ---
    if (calculatePath) {
        MatrixXi adjMatrix = calculateAdjacencyMatrix(*targetMeshPtr);

        vector<bool> isVertexInShortestPath;
        vector<bool> isEdgeInShortestPath;
        
        pair<unsigned int, double> pathResult = findShortestPathBFS(
            *targetMeshPtr,
            adjMatrix,
            startVertexId,
            endVertexId,
            isVertexInShortestPath,
            isEdgeInShortestPath
        );

        // Stampa i risultati
        if (pathResult.first > 0 || startVertexId == endVertexId) {
            cout << "\n--- Risultati Cammino Minimo (BFS) ---\n";
            cout << "Numero di lati nel cammino: " << pathResult.first << endl;
            cout << "Lunghezza totale del cammino: " << pathResult.second << endl;

            cout << "Vertici nel cammino (ID Reale): ";
            for (unsigned int i = 0; i < targetMeshPtr->Cell0DsId.size(); ++i) {
                if (isVertexInShortestPath[i]) {
                    cout << targetMeshPtr->Cell0DsId[i] << " ";
                }
            }
            cout << endl;

            cout << "Lati nel cammino (ID Reale): ";
            for (unsigned int i = 0; i < targetMeshPtr->Cell1DsId.size(); ++i) {
                if (isEdgeInShortestPath[i]) {
                    cout << targetMeshPtr->Cell1DsId[i] << " ";
                }
            }
            cout << endl;

            // Esportazione Paraview con cammino
			std::string filename = "mesh_with_path.vtk";
			PolyhedralLibrary::ExportParaview(*targetMeshPtr);		
		
        } else {
            cout << "\nNessun cammino trovato tra il vertice " << startVertexId
                 << " e il vertice " << endVertexId << ".\n";
        }
    } else {
        cout << "\nCammino minimo non richiesto (solo 5 argomenti forniti).\n";
		PolyhedralLibrary::ExportParaview(*targetMeshPtr);
    }
	

    // Scrittura su TXT
	PolyhedralLibrary::WriteCell0Ds(*targetMeshPtr);
	PolyhedralLibrary::WriteCell1Ds(*targetMeshPtr);
	PolyhedralLibrary::WriteCell2Ds(*targetMeshPtr);
	PolyhedralLibrary::WriteCell3Ds(*targetMeshPtr);
	
	*/
	
	int p = 3;
	int q = 4;
	int b = 3;
	int c = 0;
	
	PolyhedralLibrary::PolyhedralMesh mesh;
    PolyhedralLibrary::PolyhedralMesh meshTriangulated;
    PolyhedralLibrary::PolyhedralMesh meshFinal;
    PolyhedralLibrary::PolyhedralMesh meshDual;
    
	PolyhedralLibrary::generateOctahedron(mesh);
    vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
	PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
    PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
	PolyhedralLibrary::NewMesh(meshTriangulated, meshFinal, dimension);
	PolyhedralLibrary::PopulateCell3D(meshFinal, dimension);
	printMeshTriangulated(meshFinal);
	PolyhedralLibrary::CalculateDual(meshFinal, meshDual);
	PolyhedralLibrary::ProjectMeshToUnitSphere(meshDual);
	PolyhedralLibrary::ExportParaview(meshDual);

    return 0;
}

	
	