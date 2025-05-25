#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <map>
#include <algorithm> // Per std::min e std::max
#include <limits>    // Per std::numeric_limits
#include <Eigen/Dense>
#include <set>       // Per std::set per controlli efficienti di appartenenza

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
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
		vector<bool> reached(n, true);

		for (size_t i = 0; i < n; ++i) {
			if (meshTriangulated.Cell0DsFlag[i][0] == maxFlag)
				continue;
			
			for (size_t j = i+1; j < n; ++j) {
				if (meshTriangulated.Cell0DsFlag[j][0] == maxFlag)
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
						reached[i]=false;
						for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        					for (unsigned int& vertexOriginalId : meshTriangulated.Cell2DsVertices[faceId]) {
	        					if (vertexOriginalId == i){
		        					vertexOriginalId = j;
		        				}
		        			}
		        		}
						meshTriangulated.Cell0DsId[i]=j;
						break;
					}
				}
			}
			if (reached[i]) {
				meshTriangulated.Cell0DsFlag[i] = {maxFlag};
			}
		}
	}
	
	void RemoveDuplicatedEdges(PolyhedralMesh& meshTriangulated)
	{
		double tol = 1e-12;
		unsigned int maxFlag = numeric_limits<unsigned int>::max();
		size_t n = meshTriangulated.Cell1DsExtrema.rows();
		vector<bool> reached(n, true);
			
		for (size_t i = 0; i < n; ++i) {
			if (meshTriangulated.Cell1DsFlag[i] == maxFlag){
				continue; 
			}
			
			for (size_t j = i + 1; j < n; ++j) {
				if (meshTriangulated.Cell1DsFlag[j] == maxFlag){
					continue;
				}
	
				if (meshTriangulated.Cell1DsFlag[i] == meshTriangulated.Cell1DsFlag[j]) {
					int i0 = meshTriangulated.Cell1DsExtrema(i, 0);
					int i1 = meshTriangulated.Cell1DsExtrema(i, 1);
					int j0 = meshTriangulated.Cell1DsExtrema(j, 0);
					int j1 = meshTriangulated.Cell1DsExtrema(j, 1);
					
					if (((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol && (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol) || ((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol && (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol)) {
						meshTriangulated.Cell1DsFlag[i] = maxFlag;
						reached[i]=false;
						for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        					for (unsigned int& edgeOriginalId : meshTriangulated.Cell2DsEdges[faceId]) {
	        					if (edgeOriginalId == i){
		        					edgeOriginalId = j;
		        				}
		        			}
		        		}
						meshTriangulated.Cell1DsId[i]=j;
						break;
					}
				}
			}
		}	
	}
	
	/*
	
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
				// Assunzione: Cell0DsFlag[i][0] contiene l'ID del vertice unico originale a cui punta
				unsigned int unique_vertex_old_id_for_duplicate = meshTriangulated.Cell0DsFlag[i][0]; 
				
				// Reindirizza il vecchio ID duplicato al nuovo ID del suo vertice unico
				// La .at() è sicura perché unique_vertex_old_id_for_duplicate DEVE essere un ID unico già mappato
				oldToNewVertexIdMap[i] = oldToNewVertexIdMap.at(unique_vertex_old_id_for_duplicate);
			}
		}		
		
		// Ora popola le strutture finali per i vertici
		meshFinal.Cell0DsId.resize(currentNewVertexId);
		meshFinal.Cell0DsCoordinates.resize(3, currentNewVertexId);
		meshFinal.Cell0DsFlag.resize(currentNewVertexId);
	
		for (unsigned int i = 0; i < currentNewVertexId; ++i) {
			meshFinal.Cell0DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell0DsCoordinates.col(i) = temp_vertex_coordinates[i];
			meshFinal.Cell0DsFlag[i] = temp_vertex_flags[i];
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
		meshFinal.Cell1DsId.resize(currentNewEdgeId);
		meshFinal.Cell1DsExtrema.resize(currentNewEdgeId, 2);
		meshFinal.Cell1DsFlag.resize(currentNewEdgeId);
	
		for (unsigned int i = 0; i < currentNewEdgeId; ++i) {
			meshFinal.Cell1DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell1DsExtrema(i, 0) = temp_edge_extrema[i].first;
			meshFinal.Cell1DsExtrema(i, 1) = temp_edge_extrema[i].second;
			meshFinal.Cell1DsFlag[i] = temp_edge_flags[i];
		}
		
		// FACCE
		vector<vector<unsigned int>> temp_face_vertices; // Vettore temporaneo per i vertici delle facce (con i nuovi ID vertici)
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
		meshFinal.Cell2DsId.resize(currentNewFaceId);
		meshFinal.Cell2DsVertices.resize(currentNewFaceId);
		meshFinal.Cell2DsEdges.resize(currentNewFaceId);
	
		for (unsigned int i = 0; i < currentNewFaceId; ++i) {
			meshFinal.Cell2DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell2DsVertices[i] = temp_face_vertices[i];
			meshFinal.Cell2DsEdges[i] = temp_face_edges[i];
		}
		
	}
	
	void CompactPolyhedralMesh(const PolyhedralMesh& originalMesh, PolyhedralMesh& newMesh)
	{
		double tol = 1e-12;
		unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
		// -----------------------------------------------------
		// FASE 1: IDENTIFICAZIONE E MAPPATURA DEI VERTICI DUPLICATI
		// -----------------------------------------------------
		size_t num_original_vertices = originalMesh.Cell0DsCoordinates.cols();
		// Mappa old_vertex_id -> unique_old_vertex_id (ID originali)
		map<unsigned int, unsigned int> oldToUniqueVertexIdMap;
		// Vettore per tenere traccia dei vertici unici nell'ordine in cui vengono scoperti
		std::vector<unsigned int> uniqueOriginalVertexIdsInOrder;
	
		// Inizializza la mappa: ogni vertice punta a se stesso (inizialmente unico)
		for (unsigned int i = 0; i < num_original_vertices; ++i) {
			oldToUniqueVertexIdMap[i] = i;
		}
	
		for (size_t i = 0; i < num_original_vertices; ++i) {
			// Se 'i' è già stato identificato come duplicato di un vertice precedente, saltalo.
			if (oldToUniqueVertexIdMap[i] != i) {
				continue;
			}
			// 'i' è un potenziale vertice unico (il rappresentante del suo gruppo)
	
			// Prepara il set di lati per 'i' per un controllo efficiente
			std::set<unsigned int> edges_i;
			if (i < originalMesh.Cell0DsFlag.size() && !originalMesh.Cell0DsFlag[i].empty()) {
				edges_i.insert(originalMesh.Cell0DsFlag[i].begin(), originalMesh.Cell0DsFlag[i].end());
			}
	
			for (size_t j = i + 1; j < num_original_vertices; ++j) {
				// Se 'j' è già stato identificato come duplicato, saltalo.
				if (oldToUniqueVertexIdMap[j] != j) {
					continue;
				}
	
				// Controlla commonSide (se condividono almeno un lato)
				bool commonSide = false;
				if (j < originalMesh.Cell0DsFlag.size() && !originalMesh.Cell0DsFlag[j].empty()) {
					if (!edges_i.empty()) { // Solo se il vertice 'i' ha dei lati
						for (unsigned int edge_j_id : originalMesh.Cell0DsFlag[j]) {
							if (edges_i.count(edge_j_id)) {
								commonSide = true;
								break;
							}
						}
					}
				}
					
				if (commonSide) {
					// Confronto geometrico delle coordinate
					if ((originalMesh.Cell0DsCoordinates.col(i) - originalMesh.Cell0DsCoordinates.col(j)).norm() < tol) {
						// Trovato un duplicato: 'j' è geometricamente coincidente con 'i' (e su un lato comune)
						oldToUniqueVertexIdMap[j] = i; // Mappa 'j' al vertice unico 'i'
					}
				}
			}
		}
	
		// Costruisci l'elenco dei vertici unici e la loro mappatura finale (old_unique_id -> new_id)
		std::map<unsigned int, unsigned int> uniqueOldToNewVertexIdMap;
		unsigned int currentNewVertexId = 0;
		for (unsigned int i = 0; i < num_original_vertices; ++i) {
			if (oldToUniqueVertexIdMap[i] == i) { // Se 'i' è un vertice unico (mappa a se stesso)
				uniqueOldToNewVertexIdMap[i] = currentNewVertexId;
				uniqueOriginalVertexIdsInOrder.push_back(i); // Per recuperare le coordinate dopo
				currentNewVertexId++;
			}
		}
	
		// -----------------------------------------------------
		// FASE 2: POPOLAMENTO DEI VERTICI NELLA NUOVA MESH
		// -----------------------------------------------------
		newMesh.Cell0DsId.resize(currentNewVertexId);
		newMesh.Cell0DsCoordinates.resize(originalMesh.Cell0DsCoordinates.rows(), currentNewVertexId);
		newMesh.Cell0DsFlag.resize(currentNewVertexId); // Copia i flag originali dei lati
	
		for (unsigned int i = 0; i < currentNewVertexId; ++i) {
			unsigned int originalUniqueVertexId = uniqueOriginalVertexIdsInOrder[i];
			newMesh.Cell0DsId[i] = i; // Nuovi ID consecutivi
			newMesh.Cell0DsCoordinates.col(i) = originalMesh.Cell0DsCoordinates.col(originalUniqueVertexId);
			newMesh.Cell0DsFlag[i] = originalMesh.Cell0DsFlag[originalUniqueVertexId]; // Copia la lista di lati
		}
	
		// -----------------------------------------------------
		// FASE 3: IDENTIFICAZIONE E MAPPATURA DEGLI SPIGOLI DUPLICATI
		// -----------------------------------------------------
		size_t num_original_edges = originalMesh.Cell1DsExtrema.rows();
		// Mappa old_edge_id -> unique_old_edge_id (ID originali)
		std::map<unsigned int, unsigned int> oldToUniqueEdgeIdMap;
		// Vettore per tenere traccia degli spigoli unici nell'ordine in cui vengono scoperti
		std::vector<unsigned int> uniqueOriginalEdgeIdsInOrder;
	
		// Inizializza la mappa: ogni spigolo punta a se stesso (inizialmente unico)
		for (unsigned int i = 0; i < num_original_edges; ++i) {
			oldToUniqueEdgeIdMap[i] = i;
		}
	
		for (size_t i = 0; i < num_original_edges; ++i) {
			if (oldToUniqueEdgeIdMap[i] != i) { // Se 'i' è già stato marcato come duplicato, saltalo
				continue;
			}
	
			for (size_t j = i + 1; j < num_original_edges; ++j) {
				if (oldToUniqueEdgeIdMap[j] != j) { // Se 'j' è già stato marcato come duplicato, saltalo
					continue;
				}
				
				// Controlla il flag dello spigolo (se Cell1DsFlag[i] è un identificatore di "tipo" o "gruppo")
				// La tua logica originale verificava meshTriangulated.Cell1DsFlag[i] == meshTriangulated.Cell1DsFlag[j]
				// Se Cell1DsFlag non ha un significato di raggruppamento o è solo un flag generico, potresti non volerlo usare qui.
				// Per ora, lo manteniamo come nella tua funzione originale.
				if (i < originalMesh.Cell1DsFlag.size() && j < originalMesh.Cell1DsFlag.size() &&
					originalMesh.Cell1DsFlag[i] == originalMesh.Cell1DsFlag[j]) {
	
					// Recupera gli ID dei vertici degli estremi dello spigolo
					unsigned int i0_orig = originalMesh.Cell1DsExtrema(i, 0);
					unsigned int i1_orig = originalMesh.Cell1DsExtrema(i, 1);
					unsigned int j0_orig = originalMesh.Cell1DsExtrema(j, 0);
					unsigned int j1_orig = originalMesh.Cell1DsExtrema(j, 1);
	
					// IMPORTANTISSIMO: Dobbiamo confrontare le coordinate dei vertici UNICI a cui si riferiscono.
					// Usiamo oldToUniqueVertexIdMap per reindirizzare gli ID degli estremi.
					unsigned int i0_unique = oldToUniqueVertexIdMap.at(i0_orig);
					unsigned int i1_unique = oldToUniqueVertexIdMap.at(i1_orig);
					unsigned int j0_unique = oldToUniqueVertexIdMap.at(j0_orig);
					unsigned int j1_unique = oldToUniqueVertexIdMap.at(j1_orig);
	
					// Ora confrontiamo le coordinate dei vertici unici nella MESH ORIGINALE (originalMesh.Cell0DsCoordinates)
					// perché non abbiamo ancora riordinato le coordinate nella newMesh.
					// Oppure, se vuoi essere più rigoroso, potresti già aver popolato uniqueOldToNewVertexIdMap
					// e usare le coordinate da newMesh.Cell0DsCoordinates (se è già stata popolata, Fase 2).
					// Per chiarezza, uso originalMesh qui.
					Eigen::Vector3d coord_i0 = originalMesh.Cell0DsCoordinates.col(i0_unique);
					Eigen::Vector3d coord_i1 = originalMesh.Cell0DsCoordinates.col(i1_unique);
					Eigen::Vector3d coord_j0 = originalMesh.Cell0DsCoordinates.col(j0_unique);
					Eigen::Vector3d coord_j1 = originalMesh.Cell0DsCoordinates.col(j1_unique);
					
					// Verifica la coincidenza degli spigoli (due orientamenti possibili)
					if (((coord_i0 - coord_j0).norm() < tol && (coord_i1 - coord_j1).norm() < tol) ||
						((coord_i0 - coord_j1).norm() < tol && (coord_i1 - coord_j0).norm() < tol)) {
						// Trovato un duplicato: 'j' è duplicato di 'i'
						oldToUniqueEdgeIdMap[j] = i; // Mappa 'j' allo spigolo unico 'i'
					}
				}
			}
		}
	
		// Costruisci l'elenco degli spigoli unici e la loro mappatura finale (old_unique_id -> new_id)
		std::map<unsigned int, unsigned int> uniqueOldToNewEdgeIdMap;
		unsigned int currentNewEdgeId = 0;
		for (unsigned int i = 0; i < num_original_edges; ++i) {
			if (oldToUniqueEdgeIdMap[i] == i) { // Se 'i' è uno spigolo unico (mappa a se stesso)
				// Controlla se gli estremi dello spigolo unico sono degeneri dopo la compattazione dei vertici
				unsigned int v1_orig = originalMesh.Cell1DsExtrema(i, 0);
				unsigned int v2_orig = originalMesh.Cell1DsExtrema(i, 1);
				unsigned int v1_unique = oldToUniqueVertexIdMap.at(v1_orig);
				unsigned int v2_unique = oldToUniqueVertexIdMap.at(v2_orig);
	
				if (v1_unique != v2_unique) { // Se lo spigolo non è degenere (non si è ridotto a un punto)
					uniqueOldToNewEdgeIdMap[i] = currentNewEdgeId;
					uniqueOriginalEdgeIdsInOrder.push_back(i); // Per recuperare i dati dello spigolo dopo
					currentNewEdgeId++;
				} else {
					std::cerr << "Warning: Original edge " << i << " became degenerate after vertex merging and will be skipped." << std::endl;
					// Marcalo come non unico in modo che non venga usato come "padre"
					oldToUniqueEdgeIdMap[i] = maxFlag; // Un valore che indica che è stato scartato
				}
			}
		}
	
		// -----------------------------------------------------
		// FASE 4: POPOLAMENTO DEGLI SPIGOLI NELLA NUOVA MESH
		// -----------------------------------------------------
		newMesh.Cell1DsId.resize(currentNewEdgeId);
		newMesh.Cell1DsExtrema.resize(currentNewEdgeId, 2);
		newMesh.Cell1DsFlag.resize(currentNewEdgeId);
	
		for (unsigned int i = 0; i < currentNewEdgeId; ++i) {
			unsigned int originalUniqueEdgeId = uniqueOriginalEdgeIdsInOrder[i];
			newMesh.Cell1DsId[i] = i; // Nuovi ID consecutivi
			
			// Recupera gli ID dei vertici originali degli estremi dello spigolo unico
			unsigned int v1_orig = originalMesh.Cell1DsExtrema(originalUniqueEdgeId, 0);
			unsigned int v2_orig = originalMesh.Cell1DsExtrema(originalUniqueEdgeId, 1);
	
			// Mappa questi ID originali (che potrebbero essere duplicati, ma oldToUniqueVertexIdMap li reindirizza)
			// ai loro nuovi ID compattati
			newMesh.Cell1DsExtrema(i, 0) = uniqueOldToNewVertexIdMap.at(oldToUniqueVertexIdMap.at(v1_orig));
			newMesh.Cell1DsExtrema(i, 1) = uniqueOldToNewVertexIdMap.at(oldToUniqueVertexIdMap.at(v2_orig));
			
			newMesh.Cell1DsFlag[i] = originalMesh.Cell1DsFlag[originalUniqueEdgeId]; // Copia il flag dello spigolo
		}
		
		// -----------------------------------------------------
		// FASE 5: POPOLAMENTO DELLE FACCE NELLA NUOVA MESH
		// -----------------------------------------------------
		
		std::vector<unsigned int> temp_face_ids; 
		std::vector<std::vector<unsigned int>> temp_face_vertices;
		std::vector<std::vector<unsigned int>> temp_face_edges;
		
		unsigned int currentNewFaceId = 0;
	
		for (unsigned int i = 0; i < originalMesh.Cell2DsId.size(); ++i) { // Itera su tutte le facce originali
			// --- Aggiorna i vertici della faccia ---
			std::vector<unsigned int> current_face_new_vertices;
			for (unsigned int old_vertex_id : originalMesh.Cell2DsVertices[i]) {
				// Mappa l'ID del vertice originale al suo ID unico originale
				unsigned int unique_old_v_id = oldToUniqueVertexIdMap.at(old_vertex_id);
				// Mappa l'ID unico originale al suo nuovo ID compattato
				current_face_new_vertices.push_back(uniqueOldToNewVertexIdMap.at(unique_old_v_id));
			}
	
			// Verifica degenerazione vertici: meno di 3 vertici distinti
			if (current_face_new_vertices.size() < 3) {
				 std::cerr << "Warning: Skipping original Cell2D " << i << " as it has less than 3 vertices. " << std::endl;
				 continue; 
			}
			std::set<unsigned int> distinct_new_vertices(current_face_new_vertices.begin(), current_face_new_vertices.end());
			if (distinct_new_vertices.size() < 3) { 
				 std::cerr << "Warning: Skipping original Cell2D " << i << " as it became degenerate (less than 3 distinct vertices) after reindexing." << std::endl;
				 continue;
			}
	
			// --- Aggiorna gli spigoli della faccia ---
			std::vector<unsigned int> current_face_new_edges;
			bool face_has_valid_edges = true;
			if (i < originalMesh.Cell2DsEdges.size()) { // Assicurati che il vettore Cell2DsEdges esista per questo indice
				for (unsigned int old_edge_id : originalMesh.Cell2DsEdges[i]) {
					if (oldToUniqueEdgeIdMap.count(old_edge_id) > 0 && // Se l'ID originale dello spigolo è stato considerato
						oldToUniqueEdgeIdMap.at(old_edge_id) != maxFlag) { // E non è stato scartato come degenere
						
						// Mappa l'ID dello spigolo originale al suo ID unico originale
						unsigned int unique_old_edge_id = oldToUniqueEdgeIdMap.at(old_edge_id);
						// Mappa l'ID unico originale al suo nuovo ID compattato
						current_face_new_edges.push_back(uniqueOldToNewEdgeIdMap.at(unique_old_edge_id));
					} else {
						face_has_valid_edges = false; // Spigolo non valido/rimosso, la faccia è invalida
						break;
					}
				}
			} else { // Se la lista di spigoli per questa faccia non esiste
				 face_has_valid_edges = false;
			}
			
			if (!face_has_valid_edges || current_face_new_edges.empty()) { // Anche se la lista di spigoli è vuota
				std::cerr << "Warning: Skipping original Cell2D " << i << " due to referencing a removed/invalid edge or having no edges after reindexing." << std::endl;
				continue; 
			}
			
			// Se la faccia è valida dopo tutti i controlli
			temp_face_ids.push_back(currentNewFaceId); // Memorizza il nuovo ID della faccia
			temp_face_vertices.push_back(current_face_new_vertices);
			temp_face_edges.push_back(current_face_new_edges);
			currentNewFaceId++;
		}
	
		// Popola le strutture finali per le facce
		newMesh.Cell2DsId.resize(currentNewFaceId);
		newMesh.Cell2DsVertices.resize(currentNewFaceId);
		newMesh.Cell2DsEdges.resize(currentNewFaceId);
	
		for (unsigned int i = 0; i < currentNewFaceId; ++i) {
			newMesh.Cell2DsId[i] = i; // Nuovi ID consecutivi
			newMesh.Cell2DsVertices[i] = temp_face_vertices[i];
			newMesh.Cell2DsEdges[i] = temp_face_edges[i];
		}	
	}
	
	*/
}
