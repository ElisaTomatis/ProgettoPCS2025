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
	void Triangulation(int q, int b, int c, PolyhedralMesh& mesh)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;

		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
    	RemoveDuplicatedEdges(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	PopulateCell3D(meshFinal, dimension);
		ProjectMeshToUnitSphere(meshFinal);
		ExportParaview(meshFinal);
    }
    
    void TriangulationDual(int q, int b, int c, PolyhedralMesh& mesh)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;
		PolyhedralMesh meshDual;
		
		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
    	RemoveDuplicatedEdges(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	CalculateDual(meshFinal, meshDual);
		PopulateCell3D(meshDual, dimension);
		ProjectMeshToUnitSphere(meshDual);
		ExportParaview(meshDual);
    }
    
    void Triangulation2(int q, int b, int c, PolyhedralMesh& mesh)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshTriangulated2;
		PolyhedralMesh meshFinal;
		
		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
    	RemoveDuplicatedEdges(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	
    	vector<int> dimension2 = CalculateDimension2(b, q);
		triangulateAndStore2(meshFinal, meshTriangulated2, dimension2);
		PopulateCell3D(meshTriangulated2, dimension);
		ProjectMeshToUnitSphere(meshTriangulated2);
		ExportParaview(meshTriangulated2);
    }
    
    void Triangulation2Dual(int q, int b, int c, PolyhedralMesh& mesh)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshTriangulated2;
		PolyhedralMesh meshFinal;
		PolyhedralMesh meshDual;
    
		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
    	RemoveDuplicatedEdges(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	
    	vector<int> dimension2 = CalculateDimension2(b, q);
		triangulateAndStore2(meshFinal, meshTriangulated2, dimension2);
		CalculateDual(meshTriangulated2, meshDual);
		PopulateCell3D(meshDual, dimension);
		ProjectMeshToUnitSphere(meshDual);
		ExportParaview(meshDual);
    }

	void PopulateCell3D(PolyhedralMesh& meshTriangulated, const vector<int>& dimension){
		
		meshTriangulated.Cell3DsId = {0};
		meshTriangulated.NumCells0Ds = dimension[0];
		meshTriangulated.NumCells1Ds = dimension[1];
		meshTriangulated.NumCells2Ds = dimension[2];

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
