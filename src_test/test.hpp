#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Eigen"

#include <gtest/gtest.h>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace Eigen;

namespace PolyhedraTest {
	
// Dimension
TEST(TestPolyedra, TestComputePolyhedronVEF)
{
	int b = 2;
	int c = 0;
	
	// caso 1: p=3, q=3
	int q1 = 3;
	vector<int> expected1 = {10, 24, 16};
	vector<int> result1 = PolyhedralLibrary::ComputePolyhedronVEF(q1, b, c);
	EXPECT_EQ(expected1, result1);
	
	// caso 2: p=3, q=4
	int q2 = 4;
	vector<int> expected2 = {18, 48, 32};
	vector<int> result2 = PolyhedralLibrary::ComputePolyhedronVEF(q2, b, c);
	EXPECT_EQ(expected2, result2);
	
    // caso 3: p=3, q=5
	int q3 = 5;
	vector<int> expected3 = {42, 120, 80};
	vector<int> result3 = PolyhedralLibrary::ComputePolyhedronVEF(q3, b, c);
	EXPECT_EQ(expected3, result3);
}

TEST(TestPolyedra, TestCalculateDuplicated)
{
	int b = 2;
	int c = 0;
	
	// caso 1: p=3, q=3
	int q1 = 3;
	vector<int> dimension1 = {10, 24, 16};
	vector<int> expected1 = {24, 36, 16};
	vector<int> result1 = PolyhedralLibrary::CalculateDuplicated(q1, b, c, dimension1);
	EXPECT_EQ(expected1, result1);
	
	// caso 2: p=3, q=4
	int q2 = 4;
	vector<int> dimension2 = {18, 48, 32};
	vector<int> expected2 = {48, 72, 32};
	vector<int> result2 = PolyhedralLibrary::CalculateDuplicated(q2, b, c, dimension2);
	EXPECT_EQ(expected2, result2);
	
	// caso 3: p=3, q=5
	int q3 = 5;
	vector<int> dimension3 = {42, 120, 80};
	vector<int> expected3 = {120, 180 , 80};
	vector<int> result3 = PolyhedralLibrary::CalculateDuplicated(q3, b, c, dimension3);
	EXPECT_EQ(expected3, result3);
}


// Triangulation
TEST(TestPolyedra, TestTriangulationTetrahedron)
{
	PolyhedralLibrary::PolyhedralMesh meshExpected;

	// VERTICI
	meshExpected.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

	meshExpected.Cell0DsCoordinates = MatrixXd(3, 24);
	
	meshExpected.Cell0DsCoordinates.col(0)  <<  0.57735,  0.57735,  0.57735;
	meshExpected.Cell0DsCoordinates.col(1)  <<  0,        0,        0.57735;
	meshExpected.Cell0DsCoordinates.col(2)  <<  0,        0.57735,  0;
	meshExpected.Cell0DsCoordinates.col(3)  << -0.57735, -0.57735,  0.57735;
	meshExpected.Cell0DsCoordinates.col(4)  << -0.57735,  0,        0;
	meshExpected.Cell0DsCoordinates.col(5)  << -0.57735,  0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(6)  <<  0.57735,  0.57735,  0.57735;
	meshExpected.Cell0DsCoordinates.col(7)  <<  0.57735,  0,        0;
	meshExpected.Cell0DsCoordinates.col(8)  <<  0,        0,        0.57735;
	meshExpected.Cell0DsCoordinates.col(9)  <<  0.57735, -0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(10) <<  0,       -0.57735,  0;
	meshExpected.Cell0DsCoordinates.col(11) << -0.57735, -0.57735,  0.57735;
	meshExpected.Cell0DsCoordinates.col(12) << -0.57735, -0.57735,  0.57735;
	meshExpected.Cell0DsCoordinates.col(13) <<  0,       -0.57735,  0;
	meshExpected.Cell0DsCoordinates.col(14) << -0.57735,  0,        0;
	meshExpected.Cell0DsCoordinates.col(15) <<  0.57735, -0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(16) <<  0,        0,       -0.57735;
	meshExpected.Cell0DsCoordinates.col(17) << -0.57735,  0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(18) << -0.57735,  0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(19) <<  0,        0,       -0.57735;
	meshExpected.Cell0DsCoordinates.col(20) <<  0,        0.57735,  0;
	meshExpected.Cell0DsCoordinates.col(21) <<  0.57735, -0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(22) <<  0.57735,  0,        0;
	meshExpected.Cell0DsCoordinates.col(23) <<  0.57735,  0.57735,  0.57735;
	

	meshExpected.Cell0DsFlag = {
		{0, 2}, 
		{0}, 
		{2}, 
		{0, 1}, 
		{1}, 
		{1, 2}, 
		{3, 0}, 
		{3}, 
		{4294967295}, 
		{3, 4}, 
		{4}, 
		{4, 0}, 
		{4294967295}, 
		{4294967295}, 
		{4294967295},
		{4, 5}, 
		{5}, 
		{5, 1}, 
		{4294967295}, 
		{4294967295}, 
		{4294967295},
		{4294967295}, 
		{4294967295}, 
		{4294967295}
	};

	// LATI/SPIGOLI
	meshExpected.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};

	meshExpected.Cell1DsExtrema = MatrixXi(36, 2);
	meshExpected.Cell1DsExtrema << 
		0, 1, 
		1, 2, 
		2, 0,
		1, 3, 
		3, 4, 
		4, 1,
		4, 2, 
		4, 5, 
		5, 2,
		6, 7, 
		7, 8, 
		8, 6,
		7, 9, 
		9, 10, 
		10, 7,
		10, 8, 
		10, 11, 
		11, 8,
		12, 13, 
		13, 14, 
		14, 12,
		13, 15, 
		15, 16, 
		16, 13,
		16, 14, 
		16, 17, 
		17, 14,
		18, 19, 
		19, 20, 
		20, 18,
		19, 21, 
		21, 22, 
		22, 19,
		22, 20, 
		22, 23, 
		23, 20;

	meshExpected.Cell1DsFlag = {
		4294967295, 
		4294967295, 
		4294967295,
		4294967295, 
		4294967295, 
		4294967295,
		4294967295, 
		4294967295, 
		4294967295,
		4294967295, 
		4294967295, 
		0,
		4294967295, 
		4294967295, 
		4294967295,
		4294967295, 
		4294967295, 
		0,
		4, 
		4294967295, 
		1,
		4, 
		4294967295, 
		4294967295,
		4294967295, 
		4294967295, 
		1,
		5, 
		4294967295, 
		2,
		5, 
		3, 
		4294967295,
		4294967295, 
		3, 
		2
	};

	// FACCE
	meshExpected.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

	meshExpected.Cell2DsVertices = {
		{0, 1, 2}, 
		{1, 3, 4}, 
		{1, 4, 2}, 
		{2, 4, 5},
		{6, 7, 8}, 
		{7, 9, 10}, 
		{7, 10, 8}, 
		{8, 10, 11},
		{12, 13, 14}, 
		{13, 15, 16}, 
		{13, 16, 14}, 
		{14, 16, 17},
		{18, 19, 20}, 
		{19, 21, 22}, 
		{19, 22, 20}, 
		{20, 22, 23}
	};

	meshExpected.Cell2DsEdges = {
		{0, 1, 2}, 
		{3, 4, 5}, 
		{5, 6, 1}, 
		{6, 7, 8},
		{9, 10, 11}, 
		{12, 13, 14}, 
		{14, 15, 10}, 
		{15, 16, 17},
		{18, 19, 20}, 
		{21, 22, 23}, 
		{23, 24, 19}, 
		{24, 25, 26},
		{27, 28, 29}, 
		{30, 31, 32}, 
		{32, 33, 28}, 
		{33, 34, 35}
	};

	// POLIEDRI
	meshExpected.Cell3DsId = {};
	meshExpected.Cell3DsVertices = {};
	meshExpected.Cell3DsEdges = {};
	meshExpected.Cell3DsFaces = {};

	
	// mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralLibrary::PolyhedralMesh meshTriangulated;
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::generateTetrahedron(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
	PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
	PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
	
	
	
	// CONFRONTI TRA MESH TRIANGOLATA E ATTESA
	
	// Id vertici
    EXPECT_EQ(meshExpected.Cell0DsId, meshTriangulated.Cell0DsId);

    // Coordinate vertici
    EXPECT_TRUE(meshExpected.Cell0DsCoordinates.isApprox(meshTriangulated.Cell0DsCoordinates, 1e-6));

    // Flag vertici
    ASSERT_EQ(meshExpected.Cell0DsFlag.size(), meshTriangulated.Cell0DsFlag.size());
    for (size_t i = 0; i < meshExpected.Cell0DsFlag.size(); ++i) {
        EXPECT_EQ(meshExpected.Cell0DsFlag[i], meshTriangulated.Cell0DsFlag[i]) << "Mismatch at Cell0DsFlag[" << i << "]";
    }

    // Id spigoli
    EXPECT_EQ(meshExpected.Cell1DsId, meshTriangulated.Cell1DsId);

    // Estremi spigoli
    EXPECT_TRUE(meshExpected.Cell1DsExtrema == meshTriangulated.Cell1DsExtrema);

    // Flag spigoli
    EXPECT_EQ(meshExpected.Cell1DsFlag, meshTriangulated.Cell1DsFlag);

    // Id facce
    EXPECT_EQ(meshExpected.Cell2DsId, meshTriangulated.Cell2DsId);

    // Vertici delle facce
    ASSERT_EQ(meshExpected.Cell2DsVertices.size(), meshTriangulated.Cell2DsVertices.size());
    for (size_t i = 0; i < meshExpected.Cell2DsVertices.size(); ++i) {
        EXPECT_EQ(meshExpected.Cell2DsVertices[i], meshTriangulated.Cell2DsVertices[i]) << "Mismatch in Cell2DsVertices at face " << i;
    }

    // Spigoli delle facce
    ASSERT_EQ(meshExpected.Cell2DsEdges.size(), meshTriangulated.Cell2DsEdges.size());
    for (size_t i = 0; i < meshExpected.Cell2DsEdges.size(); ++i) {
        EXPECT_EQ(meshExpected.Cell2DsEdges[i], meshTriangulated.Cell2DsEdges[i]) << "Mismatch in Cell2DsEdges at face " << i;
    }

    // POLIEDRI (opzionale, se usati)
    EXPECT_EQ(meshExpected.Cell3DsId, meshTriangulated.Cell3DsId);
    EXPECT_EQ(meshExpected.Cell3DsVertices, meshTriangulated.Cell3DsVertices);
    EXPECT_EQ(meshExpected.Cell3DsEdges, meshTriangulated.Cell3DsEdges);
    EXPECT_EQ(meshExpected.Cell3DsFaces, meshTriangulated.Cell3DsFaces);
}


TEST(TestPolyedra, TestOrderedEdges)
{ 
   // mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralLibrary::PolyhedralMesh meshTriangulated;
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::generateTetrahedron(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
	PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
	PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
	
	
	// per ogni faccia i lati sono ordinati in modo che la fine dell'arco e coincida con l'inizio dell'arco successivo (e+1)%E
	// il vertice e della faccia deve corrispondere all'origine dell'arco e

    // ciclo su tutte le facce della mesh triangolata 
    for (size_t f = 0; f < meshTriangulated.Cell2DsId.size(); ++f) {
		const auto& edges = meshTriangulated.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshTriangulated.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
		
		ASSERT_EQ(vertices.size(), E) << "Numero di vertici e di lati non corrispondono per faccia " << f;
		
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			int currentEdge = edges[e]; // vertice corrente
            int nextEdge = edges[(e + 1) % E]; // vertice successivo 
			int currentEdgeOrigin = meshTriangulated.Cell1DsExtrema(currentEdge, 0);
            int currentEdgeEnd = meshTriangulated.Cell1DsExtrema(currentEdge, 1);

            int nextEdgeOrigin = meshTriangulated.Cell1DsExtrema(nextEdge, 0);

            int vertex = vertices[e];
			
			std::cout << "Face " << f << ", edge " << e << ": "
              << "currentEdgeEnd=" << currentEdgeEnd << ", "
              << "nextEdgeOrigin=" << nextEdgeOrigin << ", "
              << "vertex=" << vertex << ", "
              << "currentEdgeOrigin=" << currentEdgeOrigin << "\n";
			
			// Controllo chiusura del ciclo edge: end corrente == origin prossimo
            EXPECT_EQ(currentEdgeEnd, nextEdgeOrigin) 
                << "Edge " << currentEdge << " end non corrisponde a origin di edge " << nextEdge;
			 
            // controllo che il vertice e della faccia coincida con l'origine dell'edge corrente
            	EXPECT_EQ(vertex, currentEdgeOrigin) 
                << "Vertice " << vertex << " non coincide con origin di edge " << currentEdge;	

		 }

	}		
	
	
}

TEST(TestPolyedra, TestNotNullArea)
{
	double eps = numeric_limits<double>::epsilon();

	// mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralLibrary::PolyhedralMesh meshTriangulated;
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::generateTetrahedron(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
	PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
	PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
	
	// ciclo su tutti i triangoli
	for (size_t i = 0; i < meshTriangulated.Cell2DsVertices.size(); ++i) {
		const auto& tri = meshTriangulated.Cell2DsVertices[i];
		ASSERT_EQ(tri.size(), 3); // controllo se il triangolo ha 3 vertici
		
		// accedo alle coordinate dei vertici
		Vector3d A = meshTriangulated.Cell0DsCoordinates.col(tri[0]); // per ogni vertice prendo la colonna che contiene le coordinate
        Vector3d B = meshTriangulated.Cell0DsCoordinates.col(tri[1]);
        Vector3d C = meshTriangulated.Cell0DsCoordinates.col(tri[2]);
		
		// calcolo l'area del triangolo
		double area = 0.5 * ((B - A).cross(C - A)).norm();
		EXPECT_GT(area, eps) << "Triangolo con area nulla o quasi nulla al triangolo " << i;
	}
	
	
	
}

TEST(TestPolyedra, TestNotNullEdges){
	
	double eps = numeric_limits<double>::epsilon();

	// mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralLibrary::PolyhedralMesh meshTriangulated;
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::generateTetrahedron(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
	PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
	PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
	
	// itero su ogni lato della mesh triangolata 
	for (size_t i = 0; i < meshTriangulated.Cell1DsExtrema.rows(); ++i) {
		// prendo gli indici dei due vertici che definisco il lato i
		int vStart = meshTriangulated.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vEnd = meshTriangulated.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// prendo le coordinate dei due vertici
		Vector3d startPoint = meshTriangulated.Cell0DsCoordinates.col(vStart);
        Vector3d endPoint = meshTriangulated.Cell0DsCoordinates.col(vEnd);
		
		// calcolo la norma della lunghezza del lato
		double length = (endPoint - startPoint).norm();
		
		EXPECT_GT(length, eps) << "Lato con lunghezza nulla o quasi nulla all'edge " << i;
		
	}
}



}
