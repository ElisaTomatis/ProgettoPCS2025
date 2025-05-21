#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>

using namespace std;
namespace PolyhedralLibrary {

void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual)
{
	meshDual.Cell0DsId.resize(meshTriangulated.Cell2DsId.size());
	meshDual.Cell0DsCoordinates = MatrixXd::Zero(3, meshTriangulated.Cell2DsId.size());
	
	for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size() ; ++faceId){
		const auto& face = meshTriangulated.Cell2DsVertices[faceId];
		Vector3d barycenter = Vector3d::Zero();
		for (unsigned int vertexId : face) {
            barycenter += meshTriangulated.Cell0DsCoordinates.col(vertexId); 
        }
        barycenter /= (3.0);
        meshDual.Cell0DsCoordinates.col(faceId) = barycenter;
        meshDual.Cell0DsId[faceId]=faceId;
    }
        
        map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshTriangulated);
        // mappa che associ ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
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
		
		meshDual.Cell1DsId.resize(dualEdgesExtremaVector.size());
		meshDual.Cell1DsExtrema.resize(dualEdgesExtremaVector.size(), 2);
		for (unsigned int i = 0; i < dualEdgesExtremaVector.size(); ++i) {
			meshDual.Cell1DsId[i] = i;
			meshDual.Cell1DsExtrema(i, 0) = dualEdgesExtremaVector[i].first;
			meshDual.Cell1DsExtrema(i, 1) = dualEdgesExtremaVector[i].second;
    	}
    	 
    	map<unsigned int, vector<unsigned int>> vertexToFacesMap;
    	// mappa che per ogni vertice originale elenca le facce che lo contengono
    	for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        	for (unsigned int vertexOriginalId : meshTriangulated.Cell2DsVertices[faceId]) {
            	vertexToFacesMap[vertexOriginalId].push_back(faceId);
        	}
    	}
    	
    	vector<vector<unsigned int>> dualFacesVertices; // Vettore che conterrà i vertici per ogni nuova faccia duale
    	for (const auto& entry : vertexToFacesMap) {
        // 'entry.first' è l'ID del vertice originale.
        // 'entry.second' è il vettore degli ID delle facce originali che contengono quel vertice.
        // Questi ID di facce sono gli ID dei *vertici del duale* che formeranno la nuova faccia duale.
        	vector<unsigned int> faceVerticesForDualFace = entry.second;
        	dualFacesVertices.push_back(faceVerticesForDualFace);
        }
        
		meshDual.Cell2DsId.resize(dualFacesVertices.size());
		for(unsigned int i=0; i<dualFacesVertices.size(); ++i) {
			meshDual.Cell2DsId[i] = i;
			meshDual.Cell2DsVertices = dualFacesVertices;
		}
}

map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTrinagulated) {
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFaces;

    for (unsigned int faceId = 0; faceId < mesh.Cell2DsId.size(); ++faceId) {
        const vector<unsigned int>& faceEdges = meshTrinagulated.Cell2DsEdges[faceId];
        
        for (unsigned int edgeOriginalId : faceEdges) {
            unsigned int v1_id = meshTrinagulated.Cell1DsExtrema(edgeOriginalId, 0);
            unsigned int v2_id = meshTrinagulated.Cell1DsExtrema(edgeOriginalId, 1);
            
            pair<unsigned int, unsigned int> sortedEdgeVertices = {min(v1_id, v2_id), max(v1_id, v2_id)};
            edgeToFaces[sortedEdgeVertices].push_back(faceId);
        }
    }
    return edgeToFaces;
}
}




