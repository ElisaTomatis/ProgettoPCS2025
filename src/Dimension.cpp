#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits> 
#include <Eigen/Dense>
#include <set> 

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	void invertiValori(int& a, int& b) {
		int temp = a;
		a = b;
		b = temp;
	}
	
	vector<int> ComputePolyhedronVEF(int q, int b, int c)
	{
		vector<int> result(3); // inizializza un vettore con 3 valori, tutti -1 
	
		int T = 0;
		T = b * b + b * c + c * c;
		int V = 0;
		int E = 0;
		int F = 0;
	
		if (q == 3) {
			V = 2 * T + 2;
			E = 6 * T;
			F = 4 * T;
		}
		else if (q == 4) {
			V = 4 * T + 2;
			E = 12 * T;
			F = 8 * T;
		}
		else {
			V = 10 * T + 2;
			E = 30 * T;
			F = 20 * T;
		}
	
		result[0] = V;  // Primo elemento: V (vertici)
		result[1] = E;  // Secondo elemento: E (spigoli)
		result[2] = F;  // Terzo elemento: F (facce)
	
		return result;  // Restituisce il vettore con i valori di V, E, F
	}
	
	vector<int> CalculateDimension2(int b, const vector<int>& dimension)
	{
		vector<int> result(3);
		int V = dimension[0] + dimension[1]*(2*b-1) + dimension[2]*((3*b*b)/2 - (3*b)/2 +1);
		int E = dimension[1]*2*b + dimension[2]*((3*b*b)/2 + (3*b)/2);
		int F = dimension[2]*(3*b*b + 3*b);
		result[0] = V;  
		result[1] = E;  
		result[2] = F;
		
		return result;
	}

	vector<int> CalculateDuplicated(int q, int b, int c, const vector<int>& dimension)
	{
		vector<int> result(3);
		int subdivisionLevel = 0;
		subdivisionLevel = b + c;
		int V = dimension[0];
		int E = dimension[1];
		
		if (q == 3) {
			V += 2*4 + 6*(subdivisionLevel - 1);
			E += 6*subdivisionLevel;
		}
		else if (q == 4) {
			V += 3*6 + 12*(subdivisionLevel - 1);
			E += 12*subdivisionLevel;
		}
		else {
			V += 4*12 + 30*(subdivisionLevel - 1);
			E += 30* subdivisionLevel;
		}
		result[0] = V;  
		result[1] = E;  
		result[2] = dimension[2]; //il numero delle facce non cambia
		
		return result;
	}
	
	void RemoveDuplicatedVertices(PolyhedralMesh& meshTriangulated)
{
    double tol = 1e-12;
    unsigned int maxFlag = numeric_limits<unsigned int>::max();
    size_t n = meshTriangulated.Cell0DsCoordinates.cols();

    // Inizializza una mappa per tenere traccia del reindirizzamento finale degli ID dei vertici.
    // Inizialmente, ogni vertice punta a se stesso.
    vector<unsigned int> id_remap(n);
    for (unsigned int k = 0; k < n; ++k) {
        id_remap[k] = k;
    }

    // Un vettore per tracciare se un vertice è già stato identificato come un "master" (maxFlag)
    // o se è stato reindirizzato da un altro vertice e quindi non ha bisogno di essere processato
    // come un potenziale master.
    vector<bool> is_master_candidate(n, true);

    for (size_t i = 0; i < n; ++i) {
        // Se il vertice 'i' è già stato marcato per essere tenuto (maxFlag)
        // o se è già stato reindirizzato da un vertice precedente nel ciclo
        // (il che significa che è un duplicato e non un master), salta.
        if (meshTriangulated.Cell0DsFlag[i][0] == maxFlag || !is_master_candidate[i])
            continue;

        for (size_t j = i + 1; j < n; ++j) {
            // Se il vertice 'j' è già stato marcato per essere tenuto (maxFlag)
            // o se è già stato reindirizzato da un vertice precedente, salta.
            if (meshTriangulated.Cell0DsFlag[j][0] == maxFlag || !is_master_candidate[j])
                continue;

            bool commonSide = false;
            for (unsigned int fi : meshTriangulated.Cell0DsFlag[i]) {
                for (unsigned int fj : meshTriangulated.Cell0DsFlag[j]) {
                    if (fi == fj) {
                        commonSide = true;
                        break;
                    }
                }
                if (commonSide) break;
            }

            if (commonSide) {
                if ((meshTriangulated.Cell0DsCoordinates.col(i) - meshTriangulated.Cell0DsCoordinates.col(j)).norm() < tol) {
                    // Trovato un duplicato! Il vertice 'i' è un duplicato di 'j'.
                    // Vogliamo mantenere 'j' e reindirizzare 'i' a 'j'.

                    // Marca 'i' come NON un potenziale master.
                    is_master_candidate[i] = false;

                    // Reindirizza 'i' a 'j'.
                    // Questo è il cuore della gestione dei duplicati:
                    // tutti i vertici che precedentemente puntavano a 'i'
                    // ora devono puntare a 'j'.
                    // E 'i' stesso punterà a 'j'.
                    // Qui stiamo implementando il "Find" della Union-Find per la risoluzione.
                    unsigned int root_i = i;
                    while (id_remap[root_i] != root_i) {
                        root_i = id_remap[root_i];
                    }
                    unsigned int root_j = j;
                    while (id_remap[root_j] != root_j) {
                        root_j = id_remap[root_j];
                    }

                    // Se i "root" sono diversi, uniscili.
                    // In questo caso, stiamo dicendo che 'root_i' dovrebbe puntare a 'root_j'.
                    // Questo assicura che anche le catene di duplicati si risolvano correttamente.
                    if (root_i != root_j) {
                         id_remap[root_i] = root_j;
                    }
                    // Inoltre, assicurati che il vertice 'i' (e i suoi precedenti duplicati)
                    // puntino a 'j' o al suo master.
                    id_remap[i] = root_j;

                    // Assegna le coordinate di 'j' a 'i' (e a qualsiasi altro duplicato di 'i'
                    // che viene risolto in questa iterazione). Questo è per coerenza.
                    meshTriangulated.Cell0DsCoordinates.col(i) = meshTriangulated.Cell0DsCoordinates.col(j);
                }
            }
        }
    }

    // Fase 2: Applicare il remapping finale a tutti i vertici
    // Dobbiamo propagare le catene di reindirizzamento.
    // Un semplice loop while risolve le catene di Union-Find.
    for (unsigned int k = 0; k < n; ++k) {
        unsigned int current_id = k;
        while (id_remap[current_id] != current_id) {
            current_id = id_remap[current_id];
        }
        id_remap[k] = current_id; // Imposta il reindirizzamento finale per k
    }

    // Fase 3: Aggiornare le strutture della mesh in base al remapping finale
    for (int edgeId = 0; edgeId < meshTriangulated.Cell1DsExtrema.rows(); ++edgeId) {
		meshTriangulated.Cell1DsExtrema(edgeId, 0) = id_remap[meshTriangulated.Cell1DsExtrema(edgeId, 0)];
		meshTriangulated.Cell1DsExtrema(edgeId, 1) = id_remap[meshTriangulated.Cell1DsExtrema(edgeId, 1)];
	}
    
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int& vertexOriginalId : meshTriangulated.Cell2DsVertices[faceId]) {
            vertexOriginalId = id_remap[vertexOriginalId];
        }
    }

    // Aggiorna Cell0DsId e Cell0DsFlag in base al remapping.
    // I vertici che sono "master" avranno il loro id_remap[k] == k.
    // I vertici che sono duplicati avranno id_remap[k] != k.
    for (unsigned int k = 0; k < n; ++k) {
        meshTriangulated.Cell0DsId[k] = id_remap[k];
        if (id_remap[k] == k) {
            // Questo vertice è un master o non è stato reindirizzato
            meshTriangulated.Cell0DsFlag[k] = {maxFlag};
        } 
    }
}

	void RemoveDuplicatedEdges(PolyhedralMesh& meshTriangulated)
{
    double tol = 1e-12;
    unsigned int maxFlag = numeric_limits<unsigned int>::max();
    size_t n = meshTriangulated.Cell1DsExtrema.rows(); // Numero di lati

    // FASE 1: Inizializzazione delle strutture di reindirizzamento
    // Inizializza una mappa (vettore) per tenere traccia del reindirizzamento finale degli ID dei lati.
    // Inizialmente, ogni lato punta a se stesso.
    vector<unsigned int> id_remap(n);
    for (unsigned int k = 0; k < n; ++k) {
        id_remap[k] = k;
    }

    // Un vettore per tracciare se un lato è già stato identificato come un "master" (maxFlag)
    // o se è stato reindirizzato da un altro lato e quindi non ha bisogno di essere processato
    // come un potenziale master.
    vector<bool> is_master_candidate(n, true);

    // FASE 2: Identificazione dei duplicati e costruzione delle relazioni di reindirizzamento (Union-Find)
    for (size_t i = 0; i < n; ++i) {
        // Se il lato corrente 'i' è già stato marcato per essere tenuto (maxFlag)
        // o se è già stato reindirizzato da un lato precedente nel ciclo, salta.
        if (meshTriangulated.Cell1DsFlag[i] == maxFlag || !is_master_candidate[i])
            continue;

        for (size_t j = i + 1; j < n; ++j) {
            // Se il lato 'j' è già stato marcato per essere tenuto (maxFlag)
            // o se è già stato reindirizzato da un lato precedente, salta.
            if (meshTriangulated.Cell1DsFlag[j] == maxFlag || !is_master_candidate[j])
                continue;

            // Il tuo criterio per considerare due lati duplicati era che avessero lo stesso flag
            // E che i loro estremi (vertici) fossero gli stessi (o invertiti) entro una tolleranza.
            if (meshTriangulated.Cell1DsFlag[i] == meshTriangulated.Cell1DsFlag[j]) {
                int i0 = meshTriangulated.Cell1DsExtrema(i, 0); // Primo estremo del lato i
                int i1 = meshTriangulated.Cell1DsExtrema(i, 1); // Secondo estremo del lato i
                int j0 = meshTriangulated.Cell1DsExtrema(j, 0); // Primo estremo del lato j
                int j1 = meshTriangulated.Cell1DsExtrema(j, 1); // Secondo estremo del lato j

                // Controlla se i vertici estremi corrispondono in ordine diretto o inverso
                bool match_direct = ((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol &&
                                     (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol);

                bool match_inverse = ((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol &&
                                      (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol);

                if (match_direct || match_inverse) {
                    // Trovato un duplicato! Il lato 'i' è un duplicato del lato 'j'.
                    // Vogliamo mantenere 'j' e reindirizzare 'i' a 'j'.

                    is_master_candidate[i] = false; // Marca 'i' come NON un potenziale master.

                    // Trova il "root" (master finale) per 'i' e 'j' nella struttura id_remap
                    unsigned int root_i = i;
                    while (id_remap[root_i] != root_i) {
                        root_i = id_remap[root_i];
                    }
                    unsigned int root_j = j;
                    while (id_remap[root_j] != root_j) {
                        root_j = id_remap[root_j];
                    }

                    // Se i "root" sono diversi, uniscili.
                    // Stiamo dicendo che il master di 'i' ora deve puntare al master di 'j'.
                    // Questo garantisce che tutte le catene di duplicati si risolvano correttamente.
                    if (root_i != root_j) {
                         id_remap[root_i] = root_j; // Unisci le due componenti
                    }
                    // Assicurati che il lato 'i' stesso (e i suoi precedenti duplicati)
                    // puntino a 'j' o al suo master.
                    id_remap[i] = root_j; // Imposta il reindirizzamento diretto per 'i'

                    //meshTriangulated.Cell1DsExtrema.row(i) = meshTriangulated.Cell1DsExtrema.row(j);
                }
            }
        }
    }

    // FASE 3: Propagare le catene di reindirizzamento (Path Compression)
    // Questo ciclo assicura che ogni lato punti direttamente al suo master finale.
    for (unsigned int k = 0; k < n; ++k) {
        unsigned int current_id = k;
        while (id_remap[current_id] != current_id) {
            current_id = id_remap[current_id];
            // Optional: Path compression during traversal for future faster lookups
            id_remap[k] = current_id; // This line (if uncommented) would implement path compression during the first traversal
        }
        id_remap[k] = current_id; // Imposta il reindirizzamento finale per k
    }

    // FASE 4: Aggiornare le strutture della mesh in base al remapping finale
    // Aggiorna i riferimenti ai lati nelle facce (Cell2DsEdges)
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int& edgeOriginalId : meshTriangulated.Cell2DsEdges[faceId]) {
            edgeOriginalId = id_remap[edgeOriginalId];
        }
    }

    // Aggiorna Cell1DsId e Cell1DsFlag in base al remapping.
    for (unsigned int k = 0; k < n; ++k) {
        meshTriangulated.Cell1DsId[k] = id_remap[k]; // Ogni lato punta al suo master ID
        if (id_remap[k] == k) {
            // Questo lato è un master o non è stato reindirizzato
            meshTriangulated.Cell1DsFlag[k] = maxFlag;
        } else {
            // Questo lato è un duplicato e punterà a un master.
            // Il suo flag non viene impostato a maxFlag.
            // Potresti voler assegnare un flag specifico per i lati duplicati qui se necessario.
        }
    }
	}
	
	void NewMesh(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension)
	{
		unsigned int maxFlag = std::numeric_limits<unsigned int>::max();
		meshFinal.Cell0DsId.resize(dimension[0]);
		meshFinal.Cell0DsCoordinates = MatrixXd::Zero(3, dimension[0]);
		
		meshFinal.Cell1DsId.resize(dimension[1]);
		meshFinal.Cell1DsExtrema = MatrixXi::Zero(dimension[1], 2);
		
		meshFinal.Cell2DsId.resize(dimension[2]);
		meshFinal.Cell2DsVertices.resize(dimension[2]);
		meshFinal.Cell2DsEdges.resize(dimension[2]);
		
		unsigned int k1=0;
		unsigned int k2=0;
		
		// VERICI
    	map<unsigned int, unsigned int> oldToNewVertexIdMap; // Mappa per tradurre vecchi ID vertice -> nuovi ID vertice
    	vector<unsigned int> uniqueOldVertexIdsInOrder; // Per mantenere l'ordine dei vertici unici
    	
    	unsigned int currentNewVertexId = 0; // Contatore per i nuovi ID dei vertici
    	for (unsigned int i = 0; i < meshTriangulated.Cell0DsCoordinates.cols(); ++i) {
			if (meshTriangulated.Cell0DsFlag[i][0] == maxFlag) { 
				oldToNewVertexIdMap[i] = currentNewVertexId; // Mappa il vecchio ID unico al nuovo
				uniqueOldVertexIdsInOrder.push_back(i);     // Salva il vecchio ID per recuperare le coordinate dopo
				currentNewVertexId++; // Incrementa il contatore del nuovo ID
			}
		}
		
		// Poi, aggiungi alla mappa i vertici duplicati, reindirizzandoli ai loro vertici unici
		for (unsigned int i = 0; i < meshTriangulated.Cell0DsCoordinates.cols(); ++i) {
			// Se il vertice 'i' NON è marcato come unico (maxFlag), allora è un duplicato
			if (meshTriangulated.Cell0DsFlag[i][0] != maxFlag) {
				unsigned int uniqueVertexOldIdDuplicate = meshTriangulated.Cell0DsFlag[i][0]; 
				
				// Reindirizza il vecchio ID duplicato al nuovo ID del suo vertice unico
				// La .at() è sicura perché unique_vertex_old_id_for_duplicate DEVE essere un ID unico già mappato
				oldToNewVertexIdMap[i] = oldToNewVertexIdMap.at(uniqueVertexOldIdDuplicate);
			}
		}		
		
		// Ora popola le strutture finali per i vertici
		for (unsigned int i = 0; i < currentNewVertexId; ++i) {
			meshFinal.Cell0DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell0DsCoordinates.col(i) = meshTriangulated.Cell0DsCoordinates.col(uniqueOldVertexIdsInOrder[i]);
			meshFinal.Cell0DsFlag[i] = meshTriangulated.Cell0DsFlag[uniqueOldVertexIdsInOrder[i]];
		}
		// SPIGOLI
		map<unsigned int, unsigned int> oldToNewEdgeIdMap; // Mappa per tradurre vecchi ID spigolo -> nuovi ID spigolo
		vector<pair<unsigned int, unsigned int>> temp_edge_extrema; // Vettore temporaneo per gli estremi degli spigoli (con i nuovi ID vertici)
		vector<unsigned int> temp_edge_flags; // Vettore temporaneo per i flag degli spigoli
	
		unsigned int currentNewEdgeId = 0; // Contatore per i nuovi ID degli spigoli
		for (unsigned int i = 0; i < meshTriangulated.Cell1DsExtrema.rows(); ++i) {
			if (meshTriangulated.Cell1DsFlag[i] == maxFlag) {
				unsigned int old_v1_id = meshTriangulated.Cell1DsExtrema(i, 0);
				unsigned int old_v2_id = meshTriangulated.Cell1DsExtrema(i, 1);
	
				// Ottieni i nuovi ID dei vertici.
				// Se uno dei vertici originali non è stato mantenuto, lo spigolo non è valido.
				// oldToNewVertexIdMap.count() restituisce 1 se la chiave esiste, 0 altrimenti.
				if (oldToNewVertexIdMap.count(old_v1_id) > 0 && oldToNewVertexIdMap.count(old_v2_id) > 0) {
					unsigned int new_v1_id = oldToNewVertexIdMap[old_v1_id];
					unsigned int new_v2_id = oldToNewVertexIdMap[old_v2_id];
	
					// Ordina per consistenza (min, max) per la chiave se usata in future mappe
					temp_edge_extrema.push_back({min(new_v1_id, new_v2_id), max(new_v1_id, new_v2_id)});
					temp_edge_flags.push_back(meshTriangulated.Cell1DsFlag[i]); // Copia il flag
					oldToNewEdgeIdMap[i] = currentNewEdgeId; // Mappa il vecchio ID spigolo al nuovo
					currentNewEdgeId++;
				} else {
					// Warning se uno spigolo "valido" punta a un vertice rimosso
					cerr << "Warning: Skipping Cell1D " << i << " as one or both of its vertices were removed." << endl;
				}
			}
		}
		// Popola le strutture finali per gli spigoli

		for (unsigned int i = 0; i < currentNewEdgeId; ++i) {
			meshFinal.Cell1DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell1DsExtrema(i, 0) = temp_edge_extrema[i].first;
			meshFinal.Cell1DsExtrema(i, 1) = temp_edge_extrema[i].second;
			meshFinal.Cell1DsFlag[i] = temp_edge_flags[i];
		}
		
		vector<vector<unsigned int>> temp_face_vertices; // Vettore temporaneo per i vertici delle facce (con i nuovi ID vertici)Add commentMore actions
		vector<vector<unsigned int>> temp_face_edges; // Vettore temporaneo per gli spigoli delle facce (con i nuovi ID spigoli)
		
		unsigned int currentNewFaceId = 0; // Contatore per i nuovi ID delle facce
		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); ++i) { // Copia e aggiorna i vertici della faccia
			vector<unsigned int> current_face_new_vertices;
			bool face_valid_vertices = true;
			for (unsigned int old_vertex_id : meshTriangulated.Cell2DsVertices[i]) {
				if (oldToNewVertexIdMap.count(old_vertex_id) > 0) { // Se il vertice è stato mantenuto
					current_face_new_vertices.push_back(oldToNewVertexIdMap[old_vertex_id]);
				} else {
					face_valid_vertices = false; // Vertice non valido, la faccia è invalida
					break;
				}
			}
	
			// Copia e aggiorna gli spigoli della faccia (necessita di oldToNewEdgeIdMap)
			vector<unsigned int> current_face_new_edges;
			bool face_valid_edges = true;
			for (unsigned int old_edge_id : meshTriangulated.Cell2DsEdges[i]) {
				if (oldToNewEdgeIdMap.count(old_edge_id) > 0) { // Se lo spigolo è stato mantenuto
					current_face_new_edges.push_back(oldToNewEdgeIdMap[old_edge_id]);
				} else {
					face_valid_edges = false; // Spigolo non valido, la faccia è invalida
					break;
				}
			}
			
			// Se la faccia è valida (tutti i suoi vertici e spigoli sono stati mantenuti e riassegnati)
			if (face_valid_vertices && face_valid_edges) {
				temp_face_vertices.push_back(current_face_new_vertices);
				temp_face_edges.push_back(current_face_new_edges);
				currentNewFaceId++;
			} else {
				cerr << "Warning: Skipping Cell2D " << i << " due to invalid vertices or edges." << endl;
			}
		}
	
		// Popola le strutture finali per le facce
		for (unsigned int i = 0; i < currentNewFaceId; ++i) {
			meshFinal.Cell2DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell2DsVertices[i] = temp_face_vertices[i];
			meshFinal.Cell2DsEdges[i] = temp_face_edges[i];
		}
		
		
	}
}
