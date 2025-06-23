#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	void Triangulation(const int q, const int b, const int c, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal)
	{
		PolyhedralMesh meshTriangulated;

		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
    	RemoveDuplicatedEdges(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	PopulateCell3D(meshFinal);
    }
    
    void TriangulationDual(const int q, const int b, const int c, PolyhedralMesh& mesh, PolyhedralMesh& meshDual)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;
		
		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
    	RemoveDuplicatedEdges(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
    	CalculateDual(meshFinal, meshDual, edgeToFacesMap);
		PopulateCell3D(meshDual);
    }
    
    void Triangulation2(const int q, const int b, PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated2)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;
		
		vector<int> dimension = ComputePolyhedronVEF(q, b, 0);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, 0, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, 0,  dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
    	RemoveDuplicatedEdges(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	
    	vector<int> dimension2 = CalculateDimension2(b, q);
    	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
		triangulateAndStore2(meshFinal, meshTriangulated2, dimension2, edgeToFacesMap);
		PopulateCell3D(meshTriangulated2);
    }

	void PopulateCell3D(PolyhedralMesh& meshTriangulated){
		
		meshTriangulated.Cell3DsId = {0};
		meshTriangulated.NumCells0Ds = meshTriangulated.Cell0DsId.size();
		meshTriangulated.NumCells1Ds = meshTriangulated.Cell1DsId.size();
		meshTriangulated.NumCells2Ds = meshTriangulated.Cell2DsId.size();

		for (unsigned int i = 0; i < meshTriangulated.Cell0DsId.size(); i++){
			meshTriangulated.Cell3DsVertices.push_back(meshTriangulated.Cell0DsId[i]);
		}
		for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++){
			meshTriangulated.Cell3DsEdges.push_back(meshTriangulated.Cell1DsId[i]);

		}
		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); i++){
			meshTriangulated.Cell3DsFaces.push_back(i);
		}
	}
	

}
