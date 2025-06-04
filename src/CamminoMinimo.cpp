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


/*
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
#include <Eigen/Dense>

using namespace std;
using namespace Eigen; 
namespace PolyhedralLibrary;

// Funzione helper per calcolare la distanza euclidea tra due punti.
// Prende la mesh, e gli ID REALI dei vertici. Usa la mappa per convertirli.
double calculateDistanceById(const PolyhedralMesh& mesh, const map<unsigned int, unsigned int>& vertexIdToIndexMap, unsigned int id1, unsigned int id2) {
    
	// Converti gli ID reali in indici di colonna (0-based)
    auto it1 = vertexIdToIndexMap.find(id1);
    auto it2 = vertexIdToIndexMap.find(id2);

    if (it1 == vertexIdToIndexMap.end() || it2 == vertexIdToIndexMap.end()) {
        cerr << "Errore: ID vertice non trovato nella mappa durante il calcolo distanza.\n";
        return 0.0; // O lancia un'eccezione appropriata, e.g., std::runtime_error
    }

    unsigned int idx1 = it1->second;
    unsigned int idx2 = it2->second;

    // Accedi alle coordinate tramite gli indici di colonna
    VectorXd p1 = mesh.Cell0DsCoordinates.col(idx1);
    VectorXd p2 = mesh.Cell0DsCoordinates.col(idx2);

    return (p1 - p2).norm(); // Calcola la norma (distanza euclidea)
}

// --- Funzione per calcolare il cammino minimo usando Dijkstra ---
// startVertexId_real e endVertexId_real sono gli ID REALI dei vertici.
pair<unsigned int, double> findShortestPathDijkstra(
    PolyhedralMesh& mesh, // Passiamo la mesh per accedere e marcare i dati
    const MatrixXi& adjMatrix, // La matrice di adiacenza del grafo
    unsigned int startVertexId_real,  // ID REALE del vertice di partenza
    unsigned int endVertexId_real,    // ID REALE del vertice di arrivo
    vector<bool>& isVertexInShortestPath, // Output: vertici nel cammino (indicizzati per indice)
    vector<bool>& isEdgeInShortestPath    // Output: lati nel cammino (indicizzati per indice)
) {
    // 1. Mappa ID reale vertice -> Indice di colonna (0-based)
    map<unsigned int, unsigned int> vertexIdToIndexMap;
    for (unsigned int i = 0; i < mesh.Cell0DsId.size(); ++i) {
        vertexIdToIndexMap[mesh.Cell0DsId[i]] = i;
    }

    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols(); // Numero totale di vertici
    const unsigned int numEdgesInMesh = mesh.Cell1DsId.size();     // Numero totale di lati nella mesh

    // 2. Converti gli ID reali di partenza/arrivo in indici di colonna
    auto it_start = vertexIdToIndexMap.find(startVertexId_real);
    auto it_end = vertexIdToIndexMap.find(endVertexId_real);

    if (it_start == vertexIdToIndexMap.end() || it_end == vertexIdToIndexMap.end()) {
        cerr << "Errore: ID vertice di partenza (" << startVertexId_real
             << ") o arrivo (" << endVertexId_real << ") non valido (non trovato in Cell0DsId).\n";
        return {0, 0.0}; // Nessun cammino, 0 lati, 0 lunghezza
    }

    unsigned int startIdx = it_start->second; // Indice di colonna del vertice di partenza
    unsigned int endIdx = it_end->second;     // Indice di colonna del vertice di arrivo

    // Caso banale: partenza e arrivo sono lo stesso vertice
    if (startIdx == endIdx) {
        cout << "Il vertice di partenza e quello di arrivo sono gli stessi. Cammino minimo è 0 lati, lunghezza 0." << endl;
        return {0, 0.0};
    }

    // 3. Pre-elaborazione: Mappa per collegare una coppia di vertici (tramite INDICI)
    //    all'ID reale del lato e alla sua lunghezza. Questo serve per recuperare
    //    velocemente l'ID del lato e il suo peso durante l'esplorazione e la ricostruzione.
    map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> edgeInfoMap;
    for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
        unsigned int v1_real = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2_real = mesh.Cell1DsExtrema(i, 1);

        unsigned int v1_idx = vertexIdToIndexMap[v1_real];
        unsigned int v2_idx = vertexIdToIndexMap[v2_real];

        unsigned int edgeId = mesh.Cell1DsId[i]; // L'ID reale del lato
        double length = calculateDistanceById(mesh, vertexIdToIndexMap, v1_real, v2_real);

        // Memorizziamo l'informazione usando gli INDICI dei vertici, garantendo ordine per la chiave della mappa
        edgeInfoMap[{min(v1_idx, v2_idx), max(v1_idx, v2_idx)}] = {edgeId, length};
    }

    // 4. Inizializzazione per Dijkstra
    // dist[i]: distanza minima conosciuta dal vertice di partenza a i
    vector<double> dist(numVertices, numeric_limits<double>::infinity());
    // predVertex[i]: indice del vertice precedente nel cammino minimo a i
    vector<unsigned int> predVertex(numVertices, -1); // Usiamo -1 per indicare nessun predecessore
    // predEdge[i]: ID REALE del Cell1D usato per raggiungere i dal suo predecessore
    vector<unsigned int> predEdge(numVertices, -1);

    // visited[i]: true se il cammino più breve a 'i' è stato finalizzato
    vector<bool> visited(numVertices, false); // Nuovo: per tracciare i vertici già elaborati

    // Coda di priorità: memorizza coppie {distanza, indice_vertice}
    // std::greater per creare un min-heap (estrae l'elemento con la distanza minore)
    using QueueElem = pair<double, unsigned int>;
    priority_queue<QueueElem, vector<QueueElem>, greater<QueueElem>> pq;

    // Imposta la distanza del vertice di partenza a 0 e aggiungilo alla coda
    dist[startIdx] = 0.0;
    pq.push({0.0, startIdx});

    // 5. Esecuzione dell'algoritmo di Dijkstra
    while (!pq.empty()) {
        // Estrai il vertice 'u' con la distanza minima corrente dalla coda
        // C++17 destructuring assignment: auto [current_distance, u] = pq.top();
        // Oppure, per C++11/14:
        double current_distance = pq.top().first;
        unsigned int u = pq.top().second;
        pq.pop();

        // Se il vertice è già stato visitato, significa che abbiamo già trovato
        // il cammino più breve per esso, quindi ignoriamo questa istanza obsoleta.
        if (visited[u]) {
            continue;
        }
        visited[u] = true; // Marca il vertice come visitato/finalizzato

        // Se abbiamo raggiunto il vertice di arrivo, possiamo terminare l'algoritmo
        if (u == endIdx) {
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
                    continue; // Salta questo lato problematico
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

    // 6. Ricostruzione del cammino e calcolo delle statistiche
    // Inizializza i vettori di output con false e i marker della mesh a 0
    // isVertexInShortestPath.assign(numVertices, false);
    // isEdgeInShortestPath.assign(numEdgesInMesh, false);
	
	isVertexInShortestPath.resize(mesh.Cell0DsId.size(), false);
    isEdgeInShortestPath.resize(mesh.Cell1DsId.size(), false);
	
    mesh.Cell0DsMarker.assign(numVertices, 0); // Inizializza i marker a 0
    mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);

    unsigned int numEdgesInPath = 0;
    double totalPathLength = 0.0;

    // Se la distanza al vertice di arrivo è ancora infinito, significa che non è stato trovato alcun cammino
    if (dist[endIdx] == numeric_limits<double>::infinity()) {
        cout << "Nessun cammino trovato tra il vertice " << startVertexId_real
                  << " e il vertice " << endVertexId_real << endl;
        return {0, 0.0};
    }

    // Ricostruisci il cammino a ritroso dal vertice di arrivo al vertice di partenza
    unsigned int current_idx = endIdx; // Partiamo dall'indice del vertice di arrivo
    while (current_idx != startIdx) {
        // Marca il vertice corrente come parte del cammino minimo
        isVertexInShortestPath[current_idx] = true;
        mesh.Cell0DsMarker[current_idx] = 1;

        unsigned int prev_vertex_idx = predVertex[current_idx]; // Indice del vertice precedente nel cammino
        unsigned int edge_used_id = predEdge[current_idx];      // ID REALE del lato usato per raggiungere current_idx

        // Trova l'INDICE (posizione 0-based) del lato all'interno dei vettori Cell1DsId/Extrema
        bool foundEdgeInMesh = false;
        for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
            if (mesh.Cell1DsId[i] == edge_used_id) {
                // Marca il lato come parte del cammino minimo usando il suo INDICE nel vettore
                isEdgeInShortestPath[i] = true;
                mesh.Cell1DsMarker[i] = 1; // CORREZIONE: usa 'i' (indice), non 'edge_used_id' (ID reale)

                numEdgesInPath++;      // Incrementa il conteggio dei lati nel cammino
                foundEdgeInMesh = true;
                break;
            }
        }
        if (!foundEdgeInMesh) {
            // Questo è un errore grave: un lato dovrebbe sempre essere trovato se il predecessore esiste
            cerr << "Errore critico durante la ricostruzione: L'ID del lato (" << edge_used_id
                      << ") per il segmento (indici " << prev_vertex_idx << "->" << current_idx << ") non è stato trovato in Cell1DsId.\n";
            return {0, 0.0};
        }

        current_idx = prev_vertex_idx; // Spostati al vertice precedente e continua la ricostruzione
    }
    // Marca anche il vertice di partenza come parte del cammino
    isVertexInShortestPath[startIdx] = true;
    mesh.Cell0DsMarker[startIdx] = 1;

    // La lunghezza totale del cammino è già contenuta in dist[endIdx]
    totalPathLength = dist[endIdx];

    return {numEdgesInPath, totalPathLength};
}


// Funzione per calcolare la matrice di adiacenza (utile per BFS e Dijkstra)
MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh) {
    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols();

    // Inizializza la matrice di adiacenza con zeri
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

        // Imposta a 1 le celle corrispondenti nella matrice di adiacenza (il grafo è non diretto)
        adjMatrix(idx1, idx2) = 1;
        adjMatrix(idx2, idx1) = 1;
    }

    return adjMatrix;
}
*/


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

// Funzione helper per calcolare la distanza euclidea tra due punti.
// Prende la mesh, e gli ID REALI dei vertici. Usa la mappa per convertirli.
double calculateDistanceById(PolyhedralMesh& mesh, const map<unsigned int, unsigned int>& vertexIdToIndexMap, unsigned int id1, unsigned int id2) {
    
	// Converti gli ID reali in indici di colonna (0-based)
    auto it1 = vertexIdToIndexMap.find(id1);
    auto it2 = vertexIdToIndexMap.find(id2);

    if (it1 == vertexIdToIndexMap.end() || it2 == vertexIdToIndexMap.end()) {
        cerr << "Errore: ID vertice non trovato nella mappa durante il calcolo distanza.\n";
        return 0.0; 
    }

    unsigned int idx1 = it1->second;
    unsigned int idx2 = it2->second;

    // Accedi alle coordinate tramite gli indici di colonna
    VectorXd p1 = mesh.Cell0DsCoordinates.col(idx1);
    VectorXd p2 = mesh.Cell0DsCoordinates.col(idx2);

    return (p1 - p2).norm(); // Calcola la norma (distanza euclidea)
}


// Funzione per calcolare la matrice di adiacenza
MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh) {
    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols();

    // Inizializza la matrice di adiacenza con zeri
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

        // Imposta a 1 le celle corrispondenti nella matrice di adiacenza
        adjMatrix(idx1, idx2) = 1;
        adjMatrix(idx2, idx1) = 1;
    }

    return adjMatrix;
}

ShortestPathResult findShortestPathDijkstra(
    PolyhedralMesh& mesh,
    const MatrixXi& adjMatrix,
    unsigned int startVertexId_real,
    unsigned int endVertexId_real
) {
	const unsigned int numVertices = mesh.Cell0DsCoordinates.cols(); // Numero totale di vertici
    const unsigned int numEdgesInMesh = mesh.Cell1DsId.size();     // Numero totale di lati nella mesh

	// Inizializza il risultato, passando le dimensioni per i vettori bool
    ShortestPathResult result(0, 0.0, numVertices, numEdgesInMesh);

    // Mappa ID reale vertice -> Indice di colonna (0-based)
    map<unsigned int, unsigned int> vertexIdToIndexMap;
    for (unsigned int i = 0; i < mesh.Cell0DsId.size(); ++i) {
        vertexIdToIndexMap[mesh.Cell0DsId[i]] = i;
    }

    // Converti gli ID reali di partenza/arrivo in indici di colonna
    auto it_start = vertexIdToIndexMap.find(startVertexId_real);
    auto it_end = vertexIdToIndexMap.find(endVertexId_real);

    if (it_start == vertexIdToIndexMap.end() || it_end == vertexIdToIndexMap.end()) {
        cerr << "Errore: ID vertice di partenza (" << startVertexId_real
             << ") o arrivo (" << endVertexId_real << ") non valido (non trovato in Cell0DsId).\n";
        return result; // Nessun cammino, 0 lati, 0 lunghezza
    }

    unsigned int startIdx = it_start->second; // Indice di colonna del vertice di partenza
    unsigned int endIdx = it_end->second;     // Indice di colonna del vertice di arrivo

    // Caso banale: partenza e arrivo sono lo stesso vertice
    if (startIdx == endIdx) {
        cout << "Partenza e arrivo coincidono. Cammino nullo." << endl;
		result.verticesInPath[startIdx] = true; // Marca il vertice di partenza/arrivo
        mesh.Cell0DsMarker.assign(numVertices, 0); 
        mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
        mesh.Cell0DsMarker[startIdx] = 1; // Marca il vertice sulla mesh
        return result;
    }

    // Mappa per collegare una coppia di vertici (tramite INDICI)
    // all'ID reale del lato e alla sua lunghezza.
    map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> edgeInfoMap;
    
	for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
        unsigned int v1_real = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2_real = mesh.Cell1DsExtrema(i, 1);

        unsigned int v1_idx = vertexIdToIndexMap[v1_real];
        unsigned int v2_idx = vertexIdToIndexMap[v2_real];

        double length = calculateDistanceById(mesh, vertexIdToIndexMap, v1_real, v2_real);

        // Memorizziamo l'informazione usando gli INDICI dei vertici, garantendo ordine per la chiave della mappa
        edgeInfoMap[{min(v1_idx, v2_idx), max(v1_idx, v2_idx)}] = {mesh.Cell1DsId[i], length};
    }

    // Variabili Dijkstra
	
    // dist[i]: distanza minima conosciuta dal vertice di partenza a i
    vector<double> dist(numVertices, numeric_limits<double>::infinity());
	
    // predVertex[i]: indice del vertice precedente nel cammino minimo a i
    vector<unsigned int> predVertex(numVertices, -1); // Usiamo -1 per indicare nessun predecessore
	
    // predEdge[i]: ID REALE del Cell1D usato per raggiungere i dal suo predecessore
    vector<unsigned int> predEdge(numVertices, -1);

    // visited[i]: true se il cammino più breve a 'i' è stato finalizzato
    vector<bool> visited(numVertices, false); // Nuovo: per tracciare i vertici già elaborati

    // Coda di priorità: memorizza coppie {distanza, indice_vertice}
    // std::greater per creare un min-heap (estrae l'elemento con la distanza minore)
    using QueueElem = pair<double, unsigned int>;
    priority_queue<QueueElem, vector<QueueElem>, greater<>> pq;
	
	// priority_queue<QueueElem, vector<QueueElem>, greater<QueueElem>> pq;

    // Imposta la distanza del vertice di partenza a 0 e aggiungilo alla coda
    dist[startIdx] = 0.0;
    pq.push({0.0, startIdx});

    // algoritmo Dijkstra
    while (!pq.empty()) {
        // Estrai il vertice 'u' con la distanza minima corrente dalla coda
        // double current_distance = pq.top().first;
        // unsigned int u = pq.top().second;
		auto [current_distance, u] = pq.top();
        pq.pop();

        // Se il vertice è già stato visitato, significa che abbiamo già trovato
        // il cammino più breve per esso, quindi ignoriamo questa istanza.
        if (visited[u]) {
            continue;
        }
        visited[u] = true; // Marca il vertice come visitato/finalizzato

        // Se abbiamo raggiunto il vertice di arrivo, possiamo terminare l'algoritmo
        if (u == endIdx) {
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
    if (dist[endIdx] == numeric_limits<double>::infinity()) {
        cout << "Nessun cammino trovato tra il vertice " << startVertexId_real
                  << " e il vertice " << endVertexId_real << endl;
        return result;
    }
	
	
    // Ricostruzione del cammino e calcolo delle statistiche
    // Ricostruisci il cammino a ritroso dal vertice di arrivo al vertice di partenza
	
	// Inizializzo i marker della mesh a 0
    mesh.Cell0DsMarker.assign(numVertices, 0); // Inizializza i marker a 0
    mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
	
	result.totalLength = dist[endIdx];

    unsigned int current_idx = endIdx; // Partiamo dall'indice del vertice di arrivo
	
	while (current_idx != startIdx) {
        // Marca il vertice corrente come parte del cammino minimo
        result.verticesInPath[current_idx] = true;
        mesh.Cell0DsMarker[current_idx] = 1;

        unsigned int prev_vertex_idx = predVertex[current_idx]; // Indice del vertice precedente nel cammino
        unsigned int edge_used_id = predEdge[current_idx];      // ID REALE del lato usato per raggiungere current_idx

        // Trova l'INDICE del lato all'interno dei vettori Cell1DsId/Extrema
        bool foundEdgeInMesh = false;
        for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
            if (mesh.Cell1DsId[i] == edge_used_id) {
                // Marca il lato come parte del cammino minimo usando il suo INDICE nel vettore
                
				// mesh.Cell1DsMarker[edge_used_id] = 1; // CORREZIONE: usa 'i' (indice), non 'edge_used_id' (ID reale)
				mesh.Cell1DsMarker[i] = 1;
				
				result.edgesInPath[i] = true;
                result.numEdges++;
                
				foundEdgeInMesh = true;
                break;
            }
        }
        if (!foundEdgeInMesh) {
            // errore grave: un lato dovrebbe sempre essere trovato se il predecessore esiste
            cerr << "Errore critico durante la ricostruzione: L'ID del lato (" << edge_used_id
                      << ") per il segmento (indici " << prev_vertex_idx << "->" << current_idx << ") non è stato trovato in Cell1DsId.\n";
			// Restituisce un risultato vuoto e azzera i marker della mesh per coerenza in caso di errore grave
            mesh.Cell0DsMarker.assign(numVertices, 0);
            mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
            return ShortestPathResult(0, 0.0, numVertices, numEdgesInMesh);
        }

        current_idx = prev_vertex_idx; // Spostati al vertice precedente e continua la ricostruzione
    }
	
    // Marca anche il vertice di partenza come parte del cammino
    mesh.Cell0DsMarker[startIdx] = 1;
	result.verticesInPath[startIdx] = true;
	
	return result;
	
}

}