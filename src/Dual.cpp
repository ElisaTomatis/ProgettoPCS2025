#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>
#include <Eigen/Dense> // Assicurati di includere la libreria Eigen
#include <cmath>       // Per std::sqrt

using namespace std;
namespace PolyhedralLibrary {

void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual)
{
	// VERTICI
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
    
    // SPIGOLI
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshTriangulated);
    // mappa che associ ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
    
    unsigned int k=0;
    for (const auto& [edge_key, faces_vec] : edgeToFacesMap) {
        // Print the key (the edge)
        std::cout << "Spigolo " << k << " (" << edge_key.first << ", " << edge_key.second << ") -> Facce: [";
        
        // Print the value (the vector of faces)
        for (size_t i = 0; i < faces_vec.size(); ++i) {
            std::cout << faces_vec[i];
            if (i < faces_vec.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
        k++;
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
    
    // FACCE
    
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
    
    meshDual.Cell2DsVertices.resize(vertexToFacesMap.size());
    meshDual.Cell2DsId.resize(vertexToFacesMap.size());
    meshDual.Cell2DsEdges.resize(vertexToFacesMap.size());
    // ho tante facce duali quanti sono i vertici originali
    
    unsigned int dualFaceIdCounter = 0;
	for (const auto& entry : vertexToFacesMap) {
		unsigned int vertexOriginalId = entry.first; // Il vertice originale che è il centro della faccia duale
		std::vector<unsigned int> incidentFaces = entry.second; // Le facce originali incidenti a questo vertice
	
		// Se un vertice è incidenti a meno di 3 facce, non può formare una faccia duale valida.
		// (o se è un vertice di bordo in una mesh aperta)
		if (incidentFaces.size() < 3) {
			// Puoi decidere come gestire questi casi. Per un ottaedro chiuso, questo non dovrebbe accadere.
			// Se si tratta di una mesh aperta, questi vertici formano i bordi della faccia duale.
			// Per ora, li saltiamo o gestiamo come errore se ci aspettiamo un poliedro chiuso.
			std::cerr << "Warning: Vertice " << vertexOriginalId << " incidente a meno di 3 facce. Saltato o gestito come bordo." << std::endl;
			continue; // Salta questo vertice se non forma una faccia duale chiusa
		}
	
		// 1. Calcola il baricentro del vertice originale (per un potenziale "piano locale")
		// O semplicemente usa le coordinate del vertice stesso come origine locale
		Eigen::Vector3d vertexOriginalCoords = meshTriangulated.Cell0DsCoordinates.col(vertexOriginalId);
	
		// 2. Trova un vettore "normale" locale per il piano di proiezione (o semplicemente un vettore up arbitrario)
		// Questo è un punto delicato. Una normale valida può essere la somma normalizzata delle normali delle facce incidenti,
		// o il baricentro normalizzato dei baricentri delle facce incidenti rispetto al vertexOriginalId.
		// Per semplicità, proviamo con una normale "media" delle facce incidenti.
		Eigen::Vector3d averageNormal = Eigen::Vector3d::Zero();
		for (unsigned int faceId : incidentFaces) {
			// Per ottenere la normale della faccia, avresti bisogno di accedere ai vertici della faccia
			// e calcolare il prodotto vettoriale di due spigoli.
			// Assumiamo che tu abbia una funzione o modo per ottenere la normale di una faccia.
			// Per ora, useremo una stima meno robusta: il vettore dal vertice originale al baricentro della faccia.
			Eigen::Vector3d faceBarycenter = Eigen::Vector3d::Zero();
			const auto& faceVertices = meshTriangulated.Cell2DsVertices[faceId];
			for (unsigned int v_id : faceVertices) {
				faceBarycenter += meshTriangulated.Cell0DsCoordinates.col(v_id);
			}
			faceBarycenter /= static_cast<double>(faceVertices.size());
			averageNormal += (faceBarycenter - vertexOriginalCoords).normalized();
		}
		averageNormal.normalize(); // Normalizza la normale media
	
		// Scegli un vettore di riferimento per l'ordinamento angolare (es. un vettore arbitrario non parallelo alla normale)
		Eigen::Vector3d referenceVector;
		if (averageNormal.isApprox(Eigen::Vector3d::UnitZ())) { // Se la normale è Z, usa Y come riferimento
			referenceVector = Eigen::Vector3d::UnitY();
		} else { // Altrimenti, usa Z come riferimento (o X, non importa molto)
			referenceVector = Eigen::Vector3d::UnitZ();
		}
		// Rendi il vettore di riferimento ortogonale alla normale
		referenceVector = (referenceVector - referenceVector.dot(averageNormal) * averageNormal).normalized();
	
		// 3. Raccogli i vertici duali (baricentri delle facce originali) e i loro angoli
		std::vector<std::pair<double, unsigned int>> anglesAndDualVertices; // {angle, dual_vertex_id (original_face_id)}
	
		for (unsigned int faceId : incidentFaces) {
			// Calcola il baricentro della faccia (che è il vertice del duale)
			Eigen::Vector3d faceBarycenter = Eigen::Vector3d::Zero();
			const auto& faceVertices = meshTriangulated.Cell2DsVertices[faceId];
			for (unsigned int v_id : faceVertices) {
				faceBarycenter += meshTriangulated.Cell0DsCoordinates.col(v_id);
			}
			faceBarycenter /= static_cast<double>(faceVertices.size());
	
			// Vettore dal vertice originale al baricentro della faccia duale
			Eigen::Vector3d vecToBarycenter = faceBarycenter - vertexOriginalCoords;
	
			// Proietta il vettore sul piano normale al verticeOriginalCoords con normale `averageNormal`
			Eigen::Vector3d projectedVec = vecToBarycenter - vecToBarycenter.dot(averageNormal) * averageNormal;
			
			// Calcola l'angolo in 2D sul piano
			// Questo è il punto più critico: gli angoli devono essere robusti
			double angle = atan2(projectedVec.dot(averageNormal.cross(referenceVector)), projectedVec.dot(referenceVector));
			anglesAndDualVertices.push_back({angle, faceId});
		}
	
		// 4. Ordina i vertici duali (le facce originali) in base all'angolo
		std::sort(anglesAndDualVertices.begin(), anglesAndDualVertices.end());
	
		// 5. Costruisci la faccia duale con i vertici ordinati
		std::vector<unsigned int> orderedDualFaceVertices;
		std::vector<unsigned int> dualFaceEdges;
	
		for (size_t i = 0; i < anglesAndDualVertices.size(); ++i) {
			unsigned int currentDualVertexId = anglesAndDualVertices[i].second;
			orderedDualFaceVertices.push_back(currentDualVertexId);
	
			unsigned int nextDualVertexId;
			if (i == anglesAndDualVertices.size() - 1) {
				nextDualVertexId = anglesAndDualVertices[0].second; // Torna al primo vertice per chiudere il poligono
			} else {
				nextDualVertexId = anglesAndDualVertices[i+1].second;
			}
	
			// Trova l'ID dello spigolo duale tra currentDualVertexId e nextDualVertexId
			// Ricorda che gli ID dei vertici duali sono gli ID delle facce originali
			std::pair<unsigned int, unsigned int> dualEdgeKey = {std::min(currentDualVertexId, nextDualVertexId), std::max(currentDualVertexId, nextDualVertexId)};
			auto it_dual_edge = dualEdgeToIdMap.find(dualEdgeKey);
	
			if (it_dual_edge != dualEdgeToIdMap.end()) {
				dualFaceEdges.push_back(it_dual_edge->second);
			} else {
				std::cerr << "Errore: Spigolo duale non trovato per facce (vertici duali) "
						  << currentDualVertexId << " e " << nextDualVertexId << " attorno al vertice originale "
						  << vertexOriginalId << std::endl;
				// Questo è un errore grave: significa che la topologia non è 2-manifold
				// o che c'è un problema di precisione.
				// Potresti dover gestire un fallback o uscire.
			}
		}
		
		// Assegna la faccia duale costruita
		meshDual.Cell2DsId[dualFaceIdCounter] = dualFaceIdCounter;
		meshDual.Cell2DsVertices[dualFaceIdCounter] = orderedDualFaceVertices;
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

void ProjectMeshToUnitSphere(PolyhedralMesh& mesh) {
    for (int i = 0; i < mesh.Cell0DsCoordinates.cols(); ++i) {
        Eigen::Vector3d vertexCoords = mesh.Cell0DsCoordinates.col(i);
        double norm = vertexCoords.norm(); // Equivalente a std::sqrt(vertexCoords.squaredNorm());
        if (norm < 1e-12) { 
            cerr << "Warning: Vertice " << i << " troppo vicino all'origine. Non proiettato." << endl;
            continue; // Salta la proiezione per questo vertice
        }
        mesh.Cell0DsCoordinates.col(i) = vertexCoords / norm;
    }
}

}




