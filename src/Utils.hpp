#pragma once
#include <iostream>
#include "PolygonalMesh.hpp"
#include <Eigen/Dense>

namespace PolygonalLibrary
{
    /**
     * Triangola una faccia del poliedro e salva il risultato nella mesh triangolata.
     * 
     * mesh Mesh di input contenente le facce da triangolare.
     * meshTriangulated Mesh risultante con le facce triangolate.
     * ID della faccia da triangolare.
     * cellID ID della cella a cui appartiene la faccia.
     * dimension Informazioni di dimensione, tipicamente [numVertices, numEdges, numFaces].
     */
    void triangulateAndStore(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated,
                              unsigned int faceID, unsigned int cellID, const Eigen::Vector3i& dimension);


    /**
     * Aggiunge un lato alla mesh triangolata se non già presente.
     * 
     * a Primo estremo del lato.
     * b Secondo estremo del lato.
     * meshTriangulated Mesh in cui aggiungere il lato.
     * edgeID ID globale del nuovo lato (incrementato se un lato viene aggiunto).
     * triangleID ID del triangolo a cui il lato appartiene.
     */
    void FindAddEdge(unsigned int a, unsigned int b, PolyhedralMesh& meshTriangulated,
                     unsigned int& edgeID, unsigned int triangleID);
					 
					 
	/**
    * Esporta la mesh triangolata a paraview.
    */
	void ExportParaview(PolyhedralMesh& meshTriangulated);
}
