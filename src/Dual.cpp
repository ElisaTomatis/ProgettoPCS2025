#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>

using namespace std;
namespace PolyhedralLibrary {

void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual)
{
	// i vertici del duale sono i baricentri delle facce originali
	meshDual.Cell0DsId.resize(meshTriangulated.Cell2DsId.size());
	meshDual.Cell0DsCoordinates = MatrixXd::Zero(3, meshTriangulated.Cell2DsId.size());
	
	for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size() ; ++faceId){
		const auto& face = meshTriangulated.Cell2DsVertices[faceId]; // id dei vertici faccia originale
		Vector3d barycenter = Vector3d::Zero();
		for (unsigned int vertexId : face) {
            barycenter += meshTriangulated.Cell0DsCoordinates.col(vertexId); 
        }
        barycenter /= (3.0);
        meshDual.Cell0DsCoordinates.col(faceId) = barycenter;
        meshDual.Cell0DsId[faceId]=faceId; // l'id del vertice duale è l'id della faccia orginale
    }
        
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshTriangulated);
    // mappa che associ ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
    
    for (const auto& [edge_key, faces_vec] : edgeToFacesMap) {
        // Print the key (the edge)
        std::cout << "Spigolo (" << edge_key.first << ", " << edge_key.second << ") -> Facce: [";
        
        // Print the value (the vector of faces)
        for (size_t i = 0; i < faces_vec.size(); ++i) {
            std::cout << faces_vec[i];
            if (i < faces_vec.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }
    
    vector<pair<unsigned int, unsigned int>> dualEdgesExtremaVector;
    // Ogni spigolo interno del poliedro originale (condiviso da due facce) genera uno spigolo nel poliedro duale che connette i baricentri di quelle due facce.
    // Useremo un vettore temporaneo e poi lo convertiremo in Eigen::MatrixXi.
        
    for (const auto& entry : edgeToFacesMap) {
        const vector<unsigned int>& facesSharingEdge = entry.second; // facce
        if (facesSharingEdge.size() == 2) {
			unsigned int faceId1 = facesSharingEdge[0]; // faccia condivisa 1
			unsigned int faceId2 = facesSharingEdge[1]; // faccia condivisa 2
			pair<unsigned int, unsigned int> dualEdge = {min(faceId1, faceId2), max(faceId1, faceId2)};
			// Gli ID delle facce originali diventano gli ID dei vertici del duale
			dualEdgesExtremaVector.push_back(dualEdge);
		}
	}

	// questa struttura intermedia mi serve perchè non so a precsindere le dimensioni di meshDual.Cell1DsExtrema
	// il numero di spigoli duali è pari al numero di spigoli interni della mesh originale
		
	meshDual.Cell1DsId.resize(dualEdgesExtremaVector.size());
	meshDual.Cell1DsExtrema.resize(dualEdgesExtremaVector.size(), 2);
	for (unsigned int i = 0; i < dualEdgesExtremaVector.size(); ++i) {
		meshDual.Cell1DsId[i] = i;
		meshDual.Cell1DsExtrema(i, 0) = dualEdgesExtremaVector[i].first;
		meshDual.Cell1DsExtrema(i, 1) = dualEdgesExtremaVector[i].second;
    }
    
    map<pair<unsigned int, unsigned int>, unsigned int> dualEdgeToIdMap;
    // mappa che per ogni coppia di id di facce originali mi associa l'id dello spigolo duale
    for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); ++i) {
        unsigned int v1 = meshDual.Cell1DsExtrema(i, 0);
        unsigned int v2 = meshDual.Cell1DsExtrema(i, 1);
        dualEdgeToIdMap[{min(v1, v2), max(v1, v2)}] = i;
    }
    	 
    map<unsigned int, vector<unsigned int>> vertexToFacesMap;
    // mappa che per ogni vertice originale elenca le facce che lo contengono
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int vertexOriginalId : meshTriangulated.Cell2DsVertices[faceId]) {
            vertexToFacesMap[vertexOriginalId].push_back(faceId);
        }
    }
    
    for (const auto& [vertex_id, faces_vec] : vertexToFacesMap) {
        // Stampa la chiave (l'ID del vertice originale)
        std::cout << "Vertice " << vertex_id << " -> Facce: [";
        
        // Stampa il valore (il vettore di facce incidenti)
        for (size_t i = 0; i < faces_vec.size(); ++i) {
            std::cout << faces_vec[i];
            if (i < faces_vec.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }
    
    meshDual.Cell2DsVertices.reserve(vertexToFacesMap.size());
    meshDual.Cell2DsId.resize(vertexToFacesMap.size());
    meshDual.Cell2DsEdges.resize(vertexToFacesMap.size());
    // ho tante facce duali quanti sono i vertici originali
    
    unsigned int dualFaceIdCounter = 0;
    for (const auto& entry : vertexToFacesMap) {
        unsigned int vertexOriginalId = entry.first; // Il vertice originale che è il centro della faccia duale
        vector<unsigned int> incidentFaces = entry.second; // Le facce originali incidenti a questo vertice
                                                         // (questi saranno i vertici della faccia duale, non ordinati)

        vector<unsigned int> orderedDualFaceVertices; // Qui memorizzeremo i vertici ordinati
        vector<unsigned int> dualFaceEdges;  // Qui memorizzeremo gli spigoli della faccia duale
        vector<bool> visitedFaces(meshTriangulated.Cell2DsId.size(), false); // Per tenere traccia delle facce già visitate nel giro
        

        // Trova uno spigolo di partenza per avviare il giro
        unsigned int currentFaceId = incidentFaces[0]; // Inizia con la prima faccia incidente
        unsigned int currentEdgeIdOriginal = -1;       // ID dello spigolo originale (che ha vertexOriginalId come estremo)

        // Trova uno spigolo della currentFaceId che ha vertexOriginalId come estremo
        for (unsigned int edgeOriginalId : meshTriangulated.Cell2DsEdges[currentFaceId]) {
            unsigned int v1 = meshTriangulated.Cell1DsExtrema(edgeOriginalId, 0);
            unsigned int v2 = meshTriangulated.Cell1DsExtrema(edgeOriginalId, 1);
            if (v1 == vertexOriginalId || v2 == vertexOriginalId) {
                currentEdgeIdOriginal = edgeOriginalId;
                break;
            }
        }

        unsigned int startFaceId = currentFaceId;
        unsigned int facesVisitedCount = 0; // Contatore per il numero di facce visitate

        // Il ciclo continua finché non visitiamo tutte le facce incidenti E torniamo alla faccia iniziale
        // o finché non incontriamo una condizione di interruzione
        while (true) {
            if (visitedFaces[currentFaceId]) {
                // Abbiamo già visitato questa faccia in questo giro.
                // Se è la faccia di partenza e abbiamo visitato il numero corretto di facce, abbiamo chiuso il ciclo.
                if (currentFaceId == startFaceId && facesVisitedCount == incidentFaces.size()) {
                    break; // Ciclo chiuso correttamente
                } else {
                    // Stiamo ri-visitando una faccia ma non abbiamo chiuso il ciclo o non abbiamo visitato tutte le facce.
                    cerr << "Warning: Rilevata riconnessione inaspettata per vertice " << vertexOriginalId << endl;
                    break;
                }
            }
            visitedFaces[currentFaceId] = true;
            facesVisitedCount++;

            orderedDualFaceVertices.push_back(currentFaceId);

            const auto& facesSharingCurrentEdge = edgeToFacesMap[{min(meshTriangulated.Cell1DsExtrema(currentEdgeIdOriginal, 0), meshTriangulated.Cell1DsExtrema(currentEdgeIdOriginal, 1)), max(meshTriangulated.Cell1DsExtrema(currentEdgeIdOriginal, 0), meshTriangulated.Cell1DsExtrema(currentEdgeIdOriginal, 1))}];
            // id delle facce originali adiacenti allo spigolo currentEdgeIdOriginal
            // facesSharingCurrentEdge conterrà due ID di facce: di questi ID sarà la currentFaceId (la faccia da cui siamo "venuti"), 
            // e l'altro sarà la nextFaceId (la faccia su cui dobbiamo "andare" per continuare il giro attorno al vertexOriginalId)
            
            unsigned int nextFaceId = (unsigned int)-1;
            if (facesSharingCurrentEdge.size() == 2) {
                if (facesSharingCurrentEdge[0] == currentFaceId) {
                    nextFaceId = facesSharingCurrentEdge[1];
                } else {
                    nextFaceId = facesSharingCurrentEdge[0];
                }
            } else {
                 break; // Interrompi il giro per questa faccia duale
            }

            pair<unsigned int, unsigned int> dualEdgeKey = {min(currentFaceId, nextFaceId), max(currentFaceId, nextFaceId)};
            // chiave che identifica lo spigolo duale, lo cerchiamo nella mappa
            auto it_dual_edge = dualEdgeToIdMap.find(dualEdgeKey);
            if (it_dual_edge != dualEdgeToIdMap.end()) {
                dualFaceEdges.push_back(it_dual_edge->second); // aggiungiamo lo spigolo duale alla faccia duale
            } else {
                cerr << "Errore: Spigolo duale non trovato per facce " << currentFaceId << " e " << nextFaceId << endl;
                break;
            }

            unsigned int nextEdgeIdOriginal = (unsigned int)-1;
            // troviamo il prossimo spigolo
            for (unsigned int edgeOrigIdOnNextFace : meshTriangulated.Cell2DsEdges[nextFaceId]) {
                if (edgeOrigIdOnNextFace == currentEdgeIdOriginal) continue;
                unsigned int v1_e = meshTriangulated.Cell1DsExtrema(edgeOrigIdOnNextFace, 0);
                unsigned int v2_e = meshTriangulated.Cell1DsExtrema(edgeOrigIdOnNextFace, 1);
                if (v1_e == vertexOriginalId || v2_e == vertexOriginalId) {
	                // controlliamo se uno dei due vertici estremi dello spigolo corrente è il nostro vertexOriginalId centrale
                    nextEdgeIdOriginal = edgeOrigIdOnNextFace;
                    break;
                }
            }

            if (nextEdgeIdOriginal == (unsigned int)-1) {
                 // Questo potrebbe indicare la fine del giro se nextFaceId non ha altri spigoli
                 // incidenti a vertexOriginalId oltre a currentEdgeIdOriginal (solo per facce degeneri).
                 // In un poliedro ben formato, ci dovrebbe essere sempre un "prossimo" spigolo per continuare il giro.
                 // Se siamo arrivati qui, probabilmente il giro è completo ma non siamo tornati esattamente a startFaceId,
                 // o c'è un problema di topologia.
                 // Proviamo a fare un'ultima verifica.
                 if (nextFaceId == startFaceId && facesVisitedCount == incidentFaces.size()) {
                     // Abbiamo chiuso il ciclo correttamente con l'ultima faccia e spigolo.
                     break; 
                 }
                 cerr << "Warning: Mancante spigolo successivo per continuare il giro per vertice " << vertexOriginalId << " su faccia " << nextFaceId << endl;
                 break; 
            }
            
            currentFaceId = nextFaceId;
            currentEdgeIdOriginal = nextEdgeIdOriginal;
        }
        
        // Aggiungi la faccia duale costruita
        meshDual.Cell2DsId[dualFaceIdCounter] = dualFaceIdCounter;
        meshDual.Cell2DsVertices.push_back(orderedDualFaceVertices);

        meshDual.Cell2DsEdges[dualFaceIdCounter] = dualFaceEdges;
        dualFaceIdCounter++;
    }	
}


map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTriangulated) {
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFaces;

    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        const vector<unsigned int>& faceEdges = meshTriangulated.Cell2DsEdges[faceId]; // spigoli della faccia corrente
        
        for (unsigned int edgeOriginalId : faceEdges) {
            unsigned int v1_id = meshTriangulated.Cell1DsExtrema(edgeOriginalId, 0);
            unsigned int v2_id = meshTriangulated.Cell1DsExtrema(edgeOriginalId, 1);
            
            pair<unsigned int, unsigned int> sortedEdgeVertices = {min(v1_id, v2_id), max(v1_id, v2_id)};
            // Ordiniamo i vertici dello spigolo per avere una chiave univoca nella mappa
            edgeToFaces[sortedEdgeVertices].push_back(faceId);
        }
    }
    return edgeToFaces;
}
}




