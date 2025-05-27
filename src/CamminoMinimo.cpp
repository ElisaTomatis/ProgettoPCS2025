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

using namespace std;
using namespace Eigen;
/*
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
    const PolyhedralMesh& mesh, // Passiamo la mesh per accedere ai dati originali
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
        std::cout << "Il vertice di partenza e quello di arrivo sono gli stessi. Cammino minimo è 0 lati, lunghezza 0." << std::endl;
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
        edgeInfoMap{min(v1_idx, v2_idx), max(v1_idx, v2_idx)}] = {edgeId, length};
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
    vector<bool> isVertexInShortestPath(numVertices, false); // Vettore indicizzato per INDICE di colonna
    vector<bool> isEdgeInShortestPath(numEdgesInMesh, false); // Vettore indicizzato per INDICE in Cell1DsId

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
        unsigned int prev_vertex_idx = predecessorVertex[current_idx]; // Indice di colonna del predecessore
        unsigned int edge_used_id = predecessorEdge[current_idx];     // ID REALE del Cell1D

        // Trova l'INDICE (posizione) del Cell1D all'interno dei vettori Cell1DsId/Extrema
        bool foundEdgeIdx = false;
        for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
            if (mesh.Cell1DsId[i] == edge_used_id) {
                isEdgeInShortestPath[i] = true; // Marca l'elemento all'indice 'i'
                
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

    return {numEdgesInPath, totalPathLength};
}
}

*/