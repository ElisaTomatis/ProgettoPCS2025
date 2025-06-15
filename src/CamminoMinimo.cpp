/*
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>     // Per std::cout, std::cerr
#include <vector>       // Per std::vector
#include <queue>        // Per std::queue (per BFS)
#include <string>       // Per std::stoul
#include <cmath>        // Per std::sqrt
#include <limits>       // Per std::numeric_limits (per infinity, anche se non usato direttamente in BFS)
#include <cstdlib>      // Per EXIT_FAILURE
#include <Eigen/Dense>  // Per Eigen::MatrixXd, Eigen::VectorXd
# include <map>


using namespace std;
using namespace Eigen;
namespace PolyhedralLibrary {
	


// Funzione helper per calcolare la distanza euclidea tra due punti.
// Prende la mesh, e gli ID REALI dei vertici. Usa la mappa per convertirli.
double calculateDistanceById(
    const PolyhedralMesh& mesh,
    const map<unsigned int, unsigned int>& vertexIdToIndexMap,
    unsigned int id1,
    unsigned int id2
) {
    // Converti gli ID reali in indici di colonna
    auto it1 = vertexIdToIndexMap.find(id1);
    auto it2 = vertexIdToIndexMap.find(id2);

    if (it1 == vertexIdToIndexMap.end() || it2 == vertexIdToIndexMap.end()) {
        cerr << "Errore: ID vertice non trovato nella mappa durante il calcolo distanza.\n";
        return 0.0; // O lancia un'eccezione
    }

    unsigned int idx1 = it1->second;
    unsigned int idx2 = it2->second;

    VectorXd p1 = mesh.Cell0DsCoordinates.col(idx1);
    VectorXd p2 = mesh.Cell0DsCoordinates.col(idx2);
    return (p1 - p2).norm();
}

// --- Funzione per calcolare il cammino minimo usando BFS ---
// Ora startVertexId e endVertexId sono gli ID REALI dei vertici.
pair<unsigned int, double> findShortestPathBFS(
    PolyhedralMesh& mesh, // Passiamo la mesh per accedere ai dati originali
    const MatrixXi& adjMatrix, // La matrice di adiacenza come parametro
    unsigned int startVertexId_real,  // ID REALE del vertice di partenza
    unsigned int endVertexId_real,    // ID REALE del vertice di arrivo
    vector<bool>& isVertexInShortestPath, // Output
    vector<bool>& isEdgeInShortestPath    // Output
) {
    // 1. Costruisci la mappa ID reale -> Indice colonna
    map<unsigned int, unsigned int> vertexIdToIndexMap;
    for (unsigned int i = 0; i < mesh.Cell0DsId.size(); ++i) {
        vertexIdToIndexMap[mesh.Cell0DsId[i]] = i;
    }

    // Il numero di vertici totali (corrisponde al numero di colonne e di righe in adjMatrix)
    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols();
    const unsigned int numEdgesInMesh = mesh.Cell1DsId.size();

    // Converti gli ID reali di partenza/arrivo in indici di colonna
    auto it_start = vertexIdToIndexMap.find(startVertexId_real); // Se la chiave viene trovata, restituisce un iteratore che punta alla coppia chiave-valore trovata nella mappa.
    auto it_end = vertexIdToIndexMap.find(endVertexId_real);

	// Se la chiave non viene trovata, find() restituisce un iteratore a vertexIdToIndexMap.end()
    if (it_start == vertexIdToIndexMap.end() || it_end == vertexIdToIndexMap.end()) {
        cerr << "Errore: ID vertice di partenza o arrivo non valido (non trovato in Cell0DsId).\n";
        return {0, 0.0};
    }

    unsigned int startVertexIdx = it_start->second; // Indice di colonna del vertice di partenza
    unsigned int endVertexIdx = it_end->second;     // Indice di colonna del vertice di arrivo

    if (startVertexIdx == endVertexIdx) {
        cout << "Il vertice di partenza e quello di arrivo sono gli stessi. Cammino minimo è 0 lati, lunghezza 0." << endl;
        return {0, 0.0};
    }

    // Pre-elaborazione: Mappa per collegare una coppia di vertici (tramite INDICI) all'ID del lato e alla sua lunghezza.
    // Usiamo questa mappa per recuperare l'ID del Cell1D e la lunghezza una volta trovato il cammino.
    map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> edgeInfoMap;
    for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
        // Questi Cell1DsExtrema contengono gli ID REALI, dobbiamo convertirli in indici.
        unsigned int v1_real = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2_real = mesh.Cell1DsExtrema(i, 1);

        unsigned int v1_idx = vertexIdToIndexMap[v1_real];
        unsigned int v2_idx = vertexIdToIndexMap[v2_real];

        unsigned int edgeId = mesh.Cell1DsId[i]; // L'ID reale del lato
        double length = calculateDistanceById(mesh, vertexIdToIndexMap, v1_real, v2_real);

        // Memorizziamo l'informazione usando gli INDICI dei vertici
        edgeInfoMap[{min(v1_idx, v2_idx), max(v1_idx, v2_idx)}] = {edgeId, length};
    }

    // Inizializzazione per BFS
    queue<unsigned int> q; // La coda conterrà INDICI di vertici
    vector<bool> visited(numVertices, false);
    vector<unsigned int> predecessorVertex(numVertices, -1); // Indice del vertice precedente
    vector<unsigned int> predecessorEdge(numVertices, -1);   // ID REALE del Cell1D

    q.push(startVertexIdx);
    visited[startVertexIdx] = true;
    bool pathFound = false;

    // Esecuzione BFS
    while (!q.empty()) {
        unsigned int u_idx = q.front(); // u_idx è un INDICE di colonna
        q.pop();

        if (u_idx == endVertexIdx) {
            pathFound = true;
            break;
        }

        // Itera su tutti i possibili vicini 'w' usando la matrice di adiacenza
        for (unsigned int w_idx = 0; w_idx < numVertices; ++w_idx) { // w_idx è un INDICE di colonna
            // Se c'è un lato tra u_idx e w_idx (adjMatrix(u_idx, w_idx) == 1) E w_idx non è stato ancora visitato
            if (adjMatrix(u_idx, w_idx) == 1 && !visited[w_idx]) {
                visited[w_idx] = true;
                predecessorVertex[w_idx] = u_idx;
                
                // Trova l'ID del lato dalla mappa edgeInfoMap usando gli INDICI
                auto it = edgeInfoMap.find({min(u_idx, w_idx), max(u_idx, w_idx)});
                if (it != edgeInfoMap.end()) {
                    predecessorEdge[w_idx] = it->second.first; // L'ID REALE del Cell1D
                } else {
                    cerr << "Avviso: Lato tra indici (" << u_idx << "," << w_idx << ") presente in AdjacencyMatrix ma non trovato in edgeInfoMap.\n";
                    predecessorEdge[w_idx] = -1; // Fallback
                }
                q.push(w_idx);
            }
        }
    }

    // Ricostruzione del cammino e calcolo delle statistiche
	
	isVertexInShortestPath.assign(numVertices, false);
	isEdgeInShortestPath.assign(numEdgesInMesh, false);
    
	//vector<bool> isVertexInShortestPath(numVertices, false); // Vettore indicizzato per INDICE di colonna
    //vector<bool> isEdgeInShortestPath(numEdgesInMesh, false); // Vettore indicizzato per INDICE in Cell1DsId

	mesh.Cell0DsMarker.resize(numVertices, 0); // Inizializza i marker a 0
    mesh.Cell1DsMarker.resize(numEdgesInMesh, 0);
	
    unsigned int numEdgesInPath = 0;
    double totalPathLength = 0.0;

    if (!pathFound) {
        cout << "Nessun cammino trovato tra il vertice " << startVertexId_real
                  << " e il vertice " << endVertexId_real << endl;
        return {0, 0.0};
    }

    unsigned int current_idx = endVertexIdx; // current_idx è un INDICE di colonna
    while (current_idx != startVertexIdx) {
        isVertexInShortestPath[current_idx] = true;
		mesh.Cell0DsMarker[current_idx] = 1;
		
        unsigned int prev_vertex_idx = predecessorVertex[current_idx]; // Indice di colonna del predecessore
        unsigned int edge_used_id = predecessorEdge[current_idx];     // ID REALE del Cell1D

        // Trova l'INDICE (posizione) del Cell1D all'interno dei vettori Cell1DsId/Extrema
        bool foundEdgeIdx = false;
        for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
            if (mesh.Cell1DsId[i] == edge_used_id) {
                isEdgeInShortestPath[i] = true; // Marca l'elemento all'indice 'i'
                mesh.Cell1DsMarker[edge_used_id] = 1;
				
                // Converti gli indici di colonna in ID reali per calculateDistanceById
                unsigned int prev_vertex_real = mesh.Cell0DsId[prev_vertex_idx];
                unsigned int current_vertex_real = mesh.Cell0DsId[current_idx];

                totalPathLength += calculateDistanceById(mesh, vertexIdToIndexMap, prev_vertex_real, current_vertex_real);
                numEdgesInPath++;
                foundEdgeIdx = true;
                break;
            }
		}
        if (!foundEdgeIdx) {
            cerr << "Errore critico durante la ricostruzione: L'ID del lato (" << edge_used_id
                      << ") per il segmento (indici " << prev_vertex_idx << "->" << current_idx << ") non è stato trovato in Cell1DsId.\n";
            return {0, 0.0};
		}

        current_idx = prev_vertex_idx;
    }
    isVertexInShortestPath[startVertexIdx] = true; // Aggiungi il vertice di partenza
	mesh.Cell0DsMarker[startVertexIdx] = 1;

    return {numEdgesInPath, totalPathLength};
}


// Funzione per calcolare la matrice di adiacenza
MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh) {
    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols();
    
    // Inizializza la matrice di adiacenza con zeri
    // MatrixXi è una matrice di interi di Eigen
    MatrixXi adjMatrix = MatrixXi::Zero(numVertices, numVertices);

    // Mappa per convertire gli ID reali dei vertici in indici di colonna
    map<unsigned int, unsigned int> vertexIdToIndexMap;
    for (unsigned int i = 0; i < numVertices; ++i) {
        vertexIdToIndexMap[mesh.Cell0DsId[i]] = i;
    }

    // Itera su tutti i lati (Cell1Ds) della mesh
    for (unsigned int i = 0; i < mesh.Cell1DsId.size(); ++i) {
        // Recupera gli ID reali dei due vertici che compongono il lato
        unsigned int v1_real = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2_real = mesh.Cell1DsExtrema(i, 1);

        // Converti gli ID reali in indici di colonna
        auto it1 = vertexIdToIndexMap.find(v1_real);
        auto it2 = vertexIdToIndexMap.find(v2_real);

        if (it1 == vertexIdToIndexMap.end() || it2 == vertexIdToIndexMap.end()) {
            cerr << "Errore: ID vertice non trovato nella mappa durante la costruzione della matrice di adiacenza.\n";
            continue; // Salta questo lato e continua
        }

        unsigned int idx1 = it1->second;
        unsigned int idx2 = it2->second;

        // Imposta a 1 le celle corrispondenti nella matrice di adiacenza (il grafo non è diretto)
        adjMatrix(idx1, idx2) = 1;
        adjMatrix(idx2, idx1) = 1;
    }

    return adjMatrix;
}


}
*/



//////////////////////////////////////////////////////
// DIJKSTRA 
//////////////////////////////////////////////////////


#include "Utils.hpp"          
#include "PolyhedralMesh.hpp" 
#include <iostream>           // Per std::cout, std::cerr
#include <vector>             // Per std::vector
#include <queue>              // Per std::priority_queue
#include <map>                // Per std::map
#include <limits>             // Per std::numeric_limits
#include <cmath>              // Per std::sqrt (usato in calculateDistanceById)
#include <algorithm>          // Per std::min, std::max
#include <tuple>              // Per std::tie o destructuring assignment (C++17)
#include <utility>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen; 
namespace PolyhedralLibrary{

// IMPLEMENTAZIONE DEL COSTRUTTORE DI SHORTESTPATHRESULT
ShortestPathResult::ShortestPathResult(unsigned int nEdges, double len)
    : numEdges(nEdges), totalLength(len)
{}

// Funzione per calcolare la distanza euclidea tra due punti.
double calculateDistanceById(const PolyhedralMesh& mesh, unsigned int id1, unsigned int id2) {
    
    VectorXd p1 = mesh.Cell0DsCoordinates.col(id1);
    VectorXd p2 = mesh.Cell0DsCoordinates.col(id2);

    return (p1 - p2).norm(); // Calcola la norma (distanza euclidea)
}


// Funzione per calcolare la matrice di adiacenza
MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh) {
    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols();
	MatrixXi adjMatrix = MatrixXi::Zero(numVertices, numVertices);

    // Itera su tutti i lati (Cell1Ds) della mesh
    for (unsigned int i = 0; i < mesh.Cell1DsId.size(); ++i) {

        unsigned int v1 = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2 = mesh.Cell1DsExtrema(i, 1);
		
        adjMatrix(v1, v2) = 1;
        adjMatrix(v2, v1) = 1;
    }

    return adjMatrix;
}

ShortestPathResult findShortestPathDijkstra(
    PolyhedralMesh& mesh,
    const MatrixXi& adjMatrix,
    unsigned int startVertexId,
    unsigned int endVertexId
) {
	const unsigned int numVertices = mesh.Cell0DsCoordinates.cols(); // Numero totale di vertici
    const unsigned int numEdgesInMesh = mesh.Cell1DsId.size();     // Numero totale di lati nella mesh
    
    if (startVertexId >= numVertices) {
        cerr << "Errore: startVertexId (" << startVertexId << ") è fuori dal range valido di vertici [0, " << numVertices - 1 << "]." << endl;
        // Restituisci un risultato vuoto per indicare un errore
        return ShortestPathResult(0, 0.0);
    }
    if (endVertexId >= numVertices) {
        cerr << "Errore: endVertexId (" << endVertexId << ") è fuori dal range valido di vertici [0, " << numVertices - 1 << "]." << endl;
        // Restituisci un risultato vuoto per indicare un errore
        return ShortestPathResult(0, 0.0);
    }

	// Inizializza il risultato
    ShortestPathResult result(0, 0.0);

    // Caso banale: partenza e arrivo sono lo stesso vertice
    if (startVertexId == endVertexId) {
        cout << "Partenza e arrivo coincidono. Cammino nullo." << endl;
        mesh.Cell0DsMarker.assign(numVertices, 0); 
        mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
        mesh.Cell0DsMarker[startVertexId] = 1; // Marca il vertice sulla mesh
        return result;
    }

    // Mappa per collegare una coppia di vertici (tramite INDICI)
    // all'ID del lato e alla sua lunghezza.
    map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> edgeInfoMap;
    
	for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
        unsigned int v1 = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2 = mesh.Cell1DsExtrema(i, 1);

        double length = calculateDistanceById(mesh, v1, v2);

        // Memorizziamo l'informazione usando gli INDICI dei vertici, garantendo ordine per la chiave della mappa
        edgeInfoMap[{min(v1, v2), max(v1, v2)}] = {mesh.Cell1DsId[i], length};
    }

    // Variabili Dijkstra
	
    // dist[i]: distanza minima conosciuta dal vertice di partenza a i
    vector<double> dist(numVertices, numeric_limits<double>::infinity());
	
    // predVertex[i]: indice del vertice precedente nel cammino minimo a i
    vector<unsigned int> predVertex(numVertices, -1); // Usiamo -1 per indicare nessun predecessore
	
    // predEdge[i]: ID del Cell1D usato per raggiungere i dal suo predecessore
    vector<unsigned int> predEdge(numVertices, -1);

    // visited[i]: true se il cammino più breve a 'i' è stato finalizzato
    vector<bool> visited(numVertices, false); 

    // Coda di priorità: memorizza coppie {distanza, indice_vertice}
    // std::greater per creare un min-heap (estrae l'elemento con la distanza minore)
    using QueueElem = pair<double, unsigned int>;
    priority_queue<QueueElem, vector<QueueElem>, greater<>> pq;
	
	// priority_queue<QueueElem, vector<QueueElem>, greater<QueueElem>> pq;

    // Imposta la distanza del vertice di partenza a 0 e aggiungilo alla coda
    dist[startVertexId] = 0.0;
    pq.push({0.0, startVertexId});

    // algoritmo Dijkstra
    while (!pq.empty()) {
        // Estrai il vertice 'u' con la distanza minima corrente dalla coda
		auto [current_distance, u] = pq.top(); // accedo all'elemento
        pq.pop(); // rimuovo l'elemento

        // Se il vertice è già stato visitato, significa che abbiamo già trovato
        // il cammino più breve per esso, quindi ignoriamo questa istanza.
        if (visited[u]) {
            continue;
        }
        visited[u] = true; // Marca il vertice come visitato/finalizzato

        // Se abbiamo raggiunto il vertice di arrivo, possiamo terminare l'algoritmo
        if (u == endVertexId) {
            break;
        }

        // Itera su tutti i possibili vicini 'v' di 'u' usando la matrice di adiacenza
        for (unsigned int v = 0; v < numVertices; ++v) {
            // Se c'è un lato tra u e v (adjMatrix(u, v) == 1)
            if (adjMatrix(u, v) == 1) {
                // Recupera le informazioni sul lato (ID reale e lunghezza/peso)
                auto it_edge = edgeInfoMap.find({min(u, v), max(u, v)});
                if (it_edge == edgeInfoMap.end()) {
                    cerr << "Avviso: Lato tra indici (" << u << "," << v << ") presente in AdjacencyMatrix ma non trovato in edgeInfoMap.\n";
                    continue;
                }
                double weight = it_edge->second.second; // Peso del lato (lunghezza euclidea)

                // Operazione di "rilassamento":
                // Se la distanza calcolata a 'v' passando per 'u' è minore della distanza attuale di 'v'
                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight; // Aggiorna la distanza minima per 'v'
                    predVertex[v] = u;          // Imposta 'u' come predecessore di 'v'
                    predEdge[v] = it_edge->second.first; // Memorizza l'ID reale del lato usato
                    pq.push({dist[v], v});      // Inserisci 'v' nella coda di priorità con la nuova distanza
                }
            }
        }
    }
	
    // Se la distanza al vertice di arrivo è ancora infinito, significa che non è stato trovato alcun cammino
    if (dist[endVertexId] == numeric_limits<double>::infinity()) {
        cout << "Nessun cammino trovato tra il vertice " << startVertexId
                  << " e il vertice " << endVertexId << endl;
        mesh.Cell0DsMarker.assign(numVertices, 0);
        mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
		return result;
    }
	
	
    // Ricostruzione del cammino e calcolo delle statistiche
    // Ricostruisci il cammino a ritroso dal vertice di arrivo al vertice di partenza
	
	// Inizializzo i marker della mesh a 0
    mesh.Cell0DsMarker.assign(numVertices, 0); // Inizializza i marker a 0
    mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
	
	result.totalLength = dist[endVertexId]; // distanza totale

    unsigned int current_idx = endVertexId; // Partiamo dall'indice del vertice di arrivo
	
	while (current_idx != startVertexId) {
        // Marca il vertice corrente come parte del cammino minimo
        mesh.Cell0DsMarker[current_idx] = 1;

        unsigned int prev_vertex_idx = predVertex[current_idx]; // Indice del vertice precedente nel cammino
        unsigned int edge_used_id = predEdge[current_idx];      // ID del lato usato per raggiungere current_idx

        mesh.Cell1DsMarker[edge_used_id] = 1;
        result.numEdges++;
        
        current_idx = prev_vertex_idx; // Spostati al vertice precedente e continua la ricostruzione
    }
	
    // Marca anche il vertice di partenza come parte del cammino
    mesh.Cell0DsMarker[startVertexId] = 1;	
	
	return result;
}

}