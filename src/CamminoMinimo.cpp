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
namespace PolyhedralLibrary {
/*
// Funzione helper per calcolare la distanza euclidea tra due punti.
// Prende la mesh, e gli ID REALI dei vertici. Usa la mappa per convertirli.
double calculateDistanceById(
    const PolyhedralMesh& mesh,
    const std::map<unsigned int, unsigned int>& vertexIdToIndexMap,
    unsigned int id1,
    unsigned int id2
) {
    // Converti gli ID reali in indici di colonna
    auto it1 = vertexIdToIndexMap.find(id1);
    auto it2 = vertexIdToIndexMap.find(id2);

    if (it1 == vertexIdToIndexMap.end() || it2 == vertexIdToIndexMap.end()) {
        std::cerr << "Errore: ID vertice non trovato nella mappa durante il calcolo distanza.\n";
        return 0.0; // O lancia un'eccezione
    }

    unsigned int idx1 = it1->second;
    unsigned int idx2 = it2->second;

    Eigen::VectorXd p1 = mesh.Cell0DsCoordinates.col(idx1);
    Eigen::VectorXd p2 = mesh.Cell0DsCoordinates.col(idx2);
    return (p1 - p2).norm();
}

// --- Funzione per calcolare il cammino minimo usando BFS ---
// Ora startVertexId e endVertexId sono gli ID REALI dei vertici.
std::pair<unsigned int, double> findShortestPathBFS(
    const PolyhedralMesh& mesh, // Passiamo la mesh per accedere ai dati originali
    const Eigen::MatrixXi& adjMatrix, // La matrice di adiacenza come parametro
    unsigned int startVertexId_real,  // ID REALE del vertice di partenza
    unsigned int endVertexId_real,    // ID REALE del vertice di arrivo
    std::vector<bool>& isVertexInShortestPath, // Output
    std::vector<bool>& isEdgeInShortestPath    // Output
) {
    // 1. Costruisci la mappa ID reale -> Indice colonna
    std::map<unsigned int, unsigned int> vertexIdToIndexMap;
    for (unsigned int i = 0; i < mesh.Cell0DsId.size(); ++i) {
        vertexIdToIndexMap[mesh.Cell0DsId[i]] = i;
    }

    // Il numero di vertici totali (corrisponde al numero di colonne e di righe in adjMatrix)
    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols();
    const unsigned int numEdgesInMesh = mesh.Cell1DsId.size();

    // 0. Verifica input base
    if (numVertices == 0) {
        std::cerr << "Errore: Mesh senza vertici." << std::endl;
        return {0, 0.0};
    }
    if (adjMatrix.rows() != numVertices || adjMatrix.cols() != numVertices) {
        std::cerr << "Errore: La AdjacencyMatrix passata non corrisponde al numero di vertici della mesh.\n";
        return {0, 0.0};
    }

    // Converti gli ID reali di partenza/arrivo in indici di colonna
    auto it_start = vertexIdToIndexMap.find(startVertexId_real);
    auto it_end = vertexIdToIndexMap.find(endVertexId_real);

    if (it_start == vertexIdToIndexMap.end() || it_end == vertexIdToIndexMap.end()) {
        std::cerr << "Errore: ID vertice di partenza o arrivo non valido (non trovato in Cell0DsId).\n";
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
    std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, double>> edgeInfoMap;
    for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
        // Questi Cell1DsExtrema contengono gli ID REALI, dobbiamo convertirli in indici.
        unsigned int v1_real = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2_real = mesh.Cell1DsExtrema(i, 1);

        unsigned int v1_idx = vertexIdToIndexMap[v1_real];
        unsigned int v2_idx = vertexIdToIndexMap[v2_real];

        unsigned int edgeId = mesh.Cell1DsId[i]; // L'ID reale del lato
        double length = calculateDistanceById(mesh, vertexIdToIndexMap, v1_real, v2_real);

        // Memorizziamo l'informazione usando gli INDICI dei vertici
        edgeInfoMap[{v1_idx, v2_idx}] = {edgeId, length};
        edgeInfoMap[{v2_idx, v1_idx}] = {edgeId, length};
    }

    // Inizializzazione per BFS
    std::queue<unsigned int> q; // La coda conterrà INDICI di vertici
    std::vector<bool> visited(numVertices, false);
    std::vector<unsigned int> predecessorVertex(numVertices, -1); // Indice del vertice precedente
    std::vector<unsigned int> predecessorEdge(numVertices, -1);   // ID REALE del Cell1D

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
                auto it = edgeInfoMap.find({u_idx, w_idx});
                if (it != edgeInfoMap.end()) {
                    predecessorEdge[w_idx] = it->second.first; // L'ID REALE del Cell1D
                } else {
                    std::cerr << "Avviso: Lato tra indici (" << u_idx << "," << w_idx << ") presente in AdjacencyMatrix ma non trovato in edgeInfoMap.\n";
                    predecessorEdge[w_idx] = -1; // Fallback
                }
                q.push(w_idx);
            }
        }
    }

    // Ricostruzione del cammino e calcolo delle statistiche
    isVertexInShortestPath.assign(numVertices, false); // Vettore indicizzato per INDICE di colonna
    isEdgeInShortestPath.assign(numEdgesInMesh, false); // Vettore indicizzato per INDICE in Cell1DsId

    unsigned int numEdgesInPath = 0;
    double totalPathLength = 0.0;

    if (!pathFound) {
        std::cout << "Nessun cammino trovato tra il vertice " << startVertexId_real
                  << " e il vertice " << endVertexId_real << std::endl;
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
            std::cerr << "Errore critico durante la ricostruzione: L'ID del lato (" << edge_used_id
                      << ") per il segmento (indici " << prev_vertex_idx << "->" << current_idx << ") non è stato trovato in Cell1DsId.\n";
            return {0, 0.0};
        }

        current_idx = prev_vertex_idx;
    }
    isVertexInShortestPath[startVertexIdx] = true; // Aggiungi il vertice di partenza

    return {numEdgesInPath, totalPathLength};
}

// --- Funzione principale di esempio (main) ---
int main(int argc, char *argv[]) {
    if (argc != 5) {
        std::cerr << "Uso: " << argv[0] << " <id_vertice_1> <id_vertice_2> <b_placeholder> <c_placeholder>\n";
        return EXIT_FAILURE;
    }

    // Qui id_vertice_1 e id_vertice_2 sono gli ID REALI dei vertici che l'utente inserisce.
    unsigned int id_vertice_1_real, id_vertice_2_real;
    int b_placeholder, c_placeholder;

    try {
        id_vertice_1_real = std::stoul(argv[1]);
        id_vertice_2_real = std::stoul(argv[2]);
        b_placeholder = std::stoi(argv[3]);
        c_placeholder = std::stoi(argv[4]);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Errore: Uno o più input non sono numeri validi.\n";
        return EXIT_FAILURE;
    } catch (const std::out_of_range& e) {
        std::cerr << "Errore: Uno o più numeri sono fuori dall'intervallo valido.\n";
        return EXIT_FAILURE;
    }

    // --- Esempio di creazione di una PolyhedralMesh per test ---
    PolyhedralMesh myMesh;

    // Definiamo i vertici con ID non sequenziali, come nel tuo esempio:
    // ID reali: 21, 15, 7, 20
    // Indici di colonna in Cell0DsCoordinates: 0,  1,  2,  3
    myMesh.Cell0DsId = {21, 15, 7, 20}; // Gli ID reali nell'ordine degli indici delle colonne
    myMesh.Cell0DsCoordinates.resize(3, myMesh.Cell0DsId.size()); // 3 righe (x,y,z), N colonne

    // Assegnamo coordinate ai vertici in base al loro INDICE di colonna:
    myMesh.Cell0DsCoordinates.col(0) << 0.0, 0.0, 0.0; // Vertice ID 21 (colonna 0)
    myMesh.Cell0DsCoordinates.col(1) << 1.0, 0.0, 0.0; // Vertice ID 15 (colonna 1)
    myMesh.Cell0DsCoordinates.col(2) << 0.0, 1.0, 0.0; // Vertice ID 7  (colonna 2)
    myMesh.Cell0DsCoordinates.col(3) << 1.0, 1.0, 0.0; // Vertice ID 20 (colonna 3)


    // 5 lati (Cell1Ds)
    myMesh.Cell1DsExtrema.resize(5, 2);
    // Qui Cell1DsExtrema contiene gli ID REALI dei vertici
    myMesh.Cell1DsId.push_back(100); myMesh.Cell1DsExtrema.row(0) << 21, 15; // Lato 100: ID 21 a ID 15
    myMesh.Cell1DsId.push_back(101); myMesh.Cell1DsExtrema.row(1) << 21, 7;  // Lato 101: ID 21 a ID 7
    myMesh.Cell1DsId.push_back(102); myMesh.Cell1DsExtrema.row(2) << 15, 20; // Lato 102: ID 15 a ID 20
    myMesh.Cell1DsId.push_back(103); myMesh.Cell1DsExtrema.row(3) << 7,  20; // Lato 103: ID 7  a ID 20
    myMesh.Cell1DsId.push_back(104); myMesh.Cell1DsExtrema.row(4) << 15, 7;  // Lato 104: ID 15 a ID 7

    // --- Creazione della matrice di adiacenza ---
    // La matrice di adiacenza DEVE essere indicizzata per INDICE di colonna (0, 1, 2, 3...)
    // Non per gli ID reali (21, 15, 7, 20...).
    // Per questo, prima mappiamo gli ID reali agli indici di colonna per la costruzione.
    std::map<unsigned int, unsigned int> temp_vertexIdToIndexMap;
    for (unsigned int i = 0; i < myMesh.Cell0DsId.size(); ++i) {
        temp_vertexIdToIndexMap[myMesh.Cell0DsId[i]] = i;
    }

    Eigen::MatrixXi myAdjacencyMatrix(myMesh.Cell0DsId.size(), myMesh.Cell0DsId.size());
    myAdjacencyMatrix.setZero();

    // Popoliamo la matrice di adiacenza usando gli ID REALI da Cell1DsExtrema,
    // convertendoli in INDICI DI COLONNA.
    for (unsigned int i = 0; i < myMesh.Cell1DsExtrema.rows(); ++i) {
        unsigned int real_id_v1 = myMesh.Cell1DsExtrema(i, 0);
        unsigned int real_id_v2 = myMesh.Cell1DsExtrema(i, 1);

        unsigned int idx_v1 = temp_vertexIdToIndexMap[real_id_v1];
        unsigned int idx_v2 = temp_vertexIdToIndexMap[real_id_v2];

        myAdjacencyMatrix(idx_v1, idx_v2) = 1;
        myAdjacencyMatrix(idx_v2, idx_v1) = 1; // Grafo non diretto
    }


    std::cout << "Calcolo cammino minimo (BFS con Matrice Adiacenza) da vertice (ID reale) " << id_vertice_1_real
              << " a vertice (ID reale) " << id_vertice_2_real << std::endl;

    std::vector<bool> isVertexInPath;
    std::vector<bool> isEdgeInPath;

    std::pair<unsigned int, double> pathStats = findShortestPathBFS(
        myMesh,
        myAdjacencyMatrix,
        id_vertice_1_real, // Passiamo l'ID reale
        id_vertice_2_real, // Passiamo l'ID reale
        isVertexInPath,
        isEdgeInPath
    );

    if (pathStats.first > 0) {
        std::cout << "\nCammino minimo (BFS) trovato:\n";
        std::cout << "Numero di lati: " << pathStats.first << std::endl;
        std::cout << "Somma delle lunghezze: " << pathStats.second << std::endl;

        std::cout << "\nFLAG PER PARAVIEW (Simulazione):\n";
        std::cout << "Vertici nel cammino (ID Reale): ";
        // Stampa gli ID reali convertendo dagli indici
        for (unsigned int i = 0; i < isVertexInPath.size(); ++i) {
            if (isVertexInPath[i]) {
                std::cout << myMesh.Cell0DsId[i] << " "; // Stampa l'ID reale
            }
        }
        std::cout << std::endl;

        std::cout << "Lati nel cammino (ID Cell1D): ";
        for (unsigned int i = 0; i < isEdgeInPath.size(); ++i) {
            if (isEdgeInPath[i]) {
                std::cout << myMesh.Cell1DsId[i] << " ";
            }
        }
        std::cout << std::endl;
    }

    return 0;
}
}

*/

