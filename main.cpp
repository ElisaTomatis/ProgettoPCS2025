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
        cerr << "Uso:\n";
        cerr << "  " << argv[0] << " p q b c\n";
        cerr << "  " << argv[0] << " p q b c startVertexId endVertexId\n";
        return 1;
    }
    
    int p = stoi(argv[1]);
    int q = stoi(argv[2]);
    int b = stoi(argv[3]);
    int c = stoi(argv[4]);

	if (!((p == 3 || q == 3) && (p >= 3 && p <= 5) && (q >= 3 && q <= 5))) {
		cerr << "Errore: Almeno uno tra p e q deve essere 3. Entrambi devono essere tra 3, 4 o 5.\n";
		return 1;
	}
	
	if (!((b == c) || (b == 0) || (c == 0))) {
		cerr << "Errore: i parametri 'b' e 'c' devono essere uguali, oppure uno dei due deve essere zero.\n";
		return 1;
	}
	
	if ((b == c) && (p!=3)){
		cerr << "Errore: la triangolazione di classe II è possibile solo per i poliedri geodetici";
	}
	
    cout << "Hai inserito: p=" << p << ", q=" << q << ", b=" << b << ", c=" << c << "\n";
    if (calculatePath) {
        cout << "Cammino minimo da vertice ID: " << startVertexId << " a vertice ID: " << endVertexId << "\n";
    }
*/
    int p = 3;
    int q = 3;
    int b = 2;
    int c = 2;
	PolyhedralMesh mesh;
	PolyhedralMesh meshFinal;
    
    if (b != c){
	    
		if (p == 3 && q == 3) {
			generateTetrahedron(mesh);
			Triangulation(q, b, c, mesh, meshFinal);
			
		} else if (p == 3 && q != 3){
			if (q == 4){
				generateOctahedron(mesh);
			} else {
				generateIcosahedron(mesh);
			}
			Triangulation(q, b, c, mesh, meshFinal);
	
		} else if (q == 3 && p!= 3) {
			if (p == 4){
				generateOctahedron(mesh);
			} else {
				generateIcosahedron(mesh);
			}
			invertiValori(p, q);
			TriangulationDual(q, b, c, mesh, meshFinal);
		}
	} else {
		if (p == 3 && q == 3) {
			generateTetrahedron(mesh);
			Triangulation2(q, b, mesh, meshFinal);
			
		} else if (p == 3 && q != 3){
			if (q == 4){
				generateOctahedron(mesh);
			} else {
				generateIcosahedron(mesh);
			}
			Triangulation2(q, b, mesh, meshFinal);
		}
	} 
	/*
	if (calculatePath) {
	MatrixXi adjMatrix = calculateAdjacencyMatrix(meshFinal);
        
        ShortestPathResult pathResult = findShortestPathDijkstra(
            meshFinal,
            adjMatrix,
            startVertexId,
            endVertexId
        );

        // Stampa i risultati
        if (pathResult.numEdges > 0 || startVertexId == endVertexId) {
            cout << "\n--- Cammino Minimo ---\n";
            cout << "Numero di lati nel cammino: " << pathResult.numEdges  << endl;
            cout << "Lunghezza totale del cammino: " << pathResult.totalLength  << endl;

            // Esportazione Paraview con cammino
			ProjectMeshToUnitSphere(meshFinal);
			ExportParaview(meshFinal);		
		
        } else {
            cout << "\nNessun cammino trovato tra il vertice " << startVertexId
                 << " e il vertice " << endVertexId << ".\n";
        }
	} else {
		// Esportazione Paraview senza cammino
		ProjectMeshToUnitSphere(meshFinal);
		ExportParaview(meshFinal);
	}

	// Scrittura su TXT
	WriteCell0Ds(meshFinal);
	WriteCell1Ds(meshFinal);
	WriteCell2Ds(meshFinal);
	WriteCell3Ds(meshFinal);
	*/
    return 0;
    
}
	
	