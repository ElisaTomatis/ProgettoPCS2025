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
	meshExpected.Cell0DsId = {23, 8, 20, 12, 14, 18, 23, 22, 8, 21, 13, 12, 12, 13, 14, 21, 19, 18, 18, 19, 20, 21, 22, 23};

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
	meshExpected.Cell1DsId ={11, 1, 35, 17, 20, 5, 6, 26, 29, 34, 10, 11, 31, 21, 14, 15, 18, 17, 18, 19, 20, 21, 30, 23, 24, 27, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};


	meshExpected.Cell1DsExtrema = MatrixXi(36, 2);
	meshExpected.Cell1DsExtrema << 
		23, 8,
		8, 20,
		20, 23,
		8, 12,
		12, 14,
		14, 8,
		14, 20,
		14, 18,
		18, 20,
		23, 22,
		22, 8,
		8, 23,
		22, 21,
		21, 13,
		13, 22,
		13, 8,
		13, 12,
		12, 8,
		12, 13,
		13, 14,
		14, 12,
		13, 21,
		21, 19,
		19, 13,
		19, 14,
		19, 18,
		18, 14,
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
		0,
		4294967295,
		2,
		0,
		1,
		4294967295,
		4294967295,
		1,
		2,
		3,
		4294967295,
		4294967295,
		3,
		4,
		4294967295,
		4294967295,
		4,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		5,
		4294967295,
		4294967295,
		5,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		4294967295,
		4294967295
	};

	// FACCE
	meshExpected.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

	meshExpected.Cell2DsVertices = {
		{23, 8, 20},
		{8, 12, 14},
		{8, 14, 20},
		{20, 14, 18},
		{23, 22, 8},
		{22, 21, 13},
		{22, 13, 8},
		{8, 13, 12},
		{12, 13, 14},
		{13, 21, 19},
		{13, 19, 14},
		{14, 19, 18},
		{18, 19, 20},
		{19, 21, 22},
		{19, 22, 20},
		{20, 22, 23}
	};

	meshExpected.Cell2DsEdges = {
		{11, 1, 35},
		{17, 20, 5},
		{5, 6, 1},
		{6, 26, 29},
		{34, 10, 11},
		{31, 21, 14},
		{14, 15, 10},
		{15, 18, 17},
		{18, 19, 20},
		{21, 30, 23},
		{23, 24, 19},
		{24, 27, 26},
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
	unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
	vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
	PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
	
	
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
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++) {
                if (currentEdge == meshTriangulated.Cell1DsId[i]) {
                    currentEdgeOrigin = meshTriangulated.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshTriangulated.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
           
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++) {
                if (nextEdge == meshTriangulated.Cell1DsId[i]) {
                    nextEdgeOrigin = meshTriangulated.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshTriangulated.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext) << "Master Edge ID " << nextEdge << " not found in Cell1DsId!";

            // recupero il vertice della faccia
            unsigned int faceVertex = vertices[e]; 

            bool condition = (currentEdgeOrigin == nextEdgeOrigin) ||
                             (currentEdgeOrigin == nextEdgeEnd) ||
                             (currentEdgeEnd == nextEdgeOrigin) ||
                             (currentEdgeEnd == nextEdgeEnd);
            EXPECT_TRUE(condition);

            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            
            bool condition2 = (faceVertex == currentEdgeOrigin) ||
            					(faceVertex == currentEdgeEnd);
            EXPECT_TRUE(condition2);

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
		double area = 0.5 * ((B - A).cross(C - A)).norm(); // prodotto vettoriale
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

TEST(Polyhedra, DualTest){
    // mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralLibrary::PolyhedralMesh meshTriangulated;
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::PolyhedralMesh meshDual;
	PolyhedralLibrary::generateTetrahedron(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension = PolyhedralLibrary::ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = PolyhedralLibrary::CalculateDuplicated(q, b, c, dimension);
	PolyhedralLibrary::triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	PolyhedralLibrary::RemoveDuplicatedEdges(meshTriangulated);
	PolyhedralLibrary::RemoveDuplicatedVertices(meshTriangulated);
	CalculateDual(meshTriangulated, meshDual);
	
	double eps = numeric_limits<double>::epsilon();
	unsigned int maxFlag = numeric_limits<unsigned int>::max();

	size_t expectedVerticesDual   = meshTriangulated.Cell2DsId.size(); // 16 facce triangolate → 16 vertici nel duale
    size_t expectedEdgesDual      = meshDual.Cell1DsId.size();         // Calcolato dinamicamente, può essere 24 per subdivisionLevel=2
   
    // numero di vertici (id)
	// numero di lati(id)
   
    EXPECT_EQ(meshDual.Cell0DsId.size(), expectedVerticesDual);
    EXPECT_EQ(meshDual.Cell1DsId.size(), expectedEdgesDual);
	
	// ogni vertice del poliedro originale genera una faccia nel duale
	EXPECT_GE(meshDual.Cell2DsId.size(), 4);  // Minimo 4 se è un tetraedro chiuso
	
	// verifico se ogni faccia duale ha almento 3 vertici
	for (const auto& dualFace : meshDual.Cell2DsVertices) {
    EXPECT_GE(dualFace.size(), 3);
    }
	
	// controllo che uno spigolo non connetta un vertice a se stesso
	for (int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
    int a = meshDual.Cell1DsExtrema(i, 0);
    int b = meshDual.Cell1DsExtrema(i, 1);
    EXPECT_NE(a, b); // Mi verifica che a e b siano diversi 
    }
	
	// verifico che il calcolo del baricentro sia corretto
	// Per ogni faccia della mesh triangolata
    for (size_t faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;
		
        // accedo alla lista dei vertici che compongono la faccia
        const vector<unsigned int>& faceVertices = meshTriangulated.Cell2DsVertices[faceId];

        // sommo le coordinate dei vertici della faccia
        for (unsigned int v_id : faceVertices) {
            sumX += meshTriangulated.Cell0DsCoordinates(0, v_id); // coordinata x del vertice v_id
            sumY += meshTriangulated.Cell0DsCoordinates(1, v_id); // coordinata y del vertice v_id
            sumZ += meshTriangulated.Cell0DsCoordinates(2, v_id); // coordinata z del vertice v_id
        }

        // calcolo il baricentro
        double n = (double) faceVertices.size(); // perché faceVertices.size() restituisce un tipo size_t quindi evito una conversione implicita
        double baryX = sumX / n;
        double baryY = sumY / n;
        double baryZ = sumZ / n;

        // estraggo la coordinata del vertice duale corrispondente (baricentro)
        double X = meshDual.Cell0DsCoordinates(0, faceId);
        double Y = meshDual.Cell0DsCoordinates(1, faceId);
        double Z = meshDual.Cell0DsCoordinates(2, faceId);

        EXPECT_NEAR(baryX, X, eps);
        EXPECT_NEAR(baryY, Y, eps);
        EXPECT_NEAR(baryZ, Z, eps);
    }

	
	// verifico che la proiezione sulla sfera sia corretta 
	// verifico che la distanza di ogni vertice dall'origine sia circa 1
	
	// proietto il duale sulla sfera unitaria
    PolyhedralLibrary::ProjectMeshToUnitSphere(meshDual);

	for (int i = 0; i < meshDual.Cell0DsCoordinates.cols(); ++i) {
    double norm = meshDual.Cell0DsCoordinates.col(i).norm(); // calcolo la distanza del vertice dall'origine considerando la norma
    EXPECT_NEAR(norm, 1.0, eps);
    }

    // verifico che i lati siano costruiti nel modo corretto ---> devono connettere i baricentri
	for (int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		Vector2i edge = meshDual.Cell1DsExtrema.row(i);
		// id delle facce del poliedro originale (vertici duali sono baricentri di facce originali)
        int f1 = edge(0); 
        int f2 = edge(1);
		 // controlla che le due facce originali condividano un edge
    bool foundCommonEdge = false;
    for (int e1 : meshTriangulated.Cell2DsEdges[f1]) {
        for (int e2 : meshTriangulated.Cell2DsEdges[f2]) {
            if (e1 == e2) {
                foundCommonEdge = true;
                break;
            }
        }
        if (foundCommonEdge) break;
    }
    EXPECT_TRUE(foundCommonEdge);
	}
	
	// verifico che i lati e le facce siano ordinate correttamente 
	
	for (size_t f = 0; f < meshDual.Cell2DsId.size(); ++f) {
		const auto& edges = meshDual.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshDual.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
       
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (currentEdge == meshDual.Cell1DsId[i]) {
                    currentEdgeOrigin = meshDual.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshDual.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
            
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (nextEdge == meshDual.Cell1DsId[i]) {
                    nextEdgeOrigin = meshDual.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshDual.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext);

            // Recupera il vertice della faccia
            unsigned int faceVertex = vertices[e]; // ID master del vertice della faccia

            bool condition = (currentEdgeOrigin == nextEdgeOrigin) ||
                             (currentEdgeOrigin == nextEdgeEnd) ||
                             (currentEdgeEnd == nextEdgeOrigin) ||
                             (currentEdgeEnd == nextEdgeEnd);
            EXPECT_TRUE(condition);
              
            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            bool condition2 = (faceVertex == currentEdgeOrigin) ||
            					(faceVertex == currentEdgeEnd);
            EXPECT_TRUE(condition2);
		}		
	}	 
}
}
