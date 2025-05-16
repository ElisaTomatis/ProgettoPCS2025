#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "UCDUtilities.hpp"
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

namespace PolyhedralLibrary
{
	void ExportParaview(const PolyhedralMesh& meshTriangulated){
		
		Gedim::UCDUtilities utilities;
		{
			utilities.ExportPoints("./Cell0Ds.inp",
									 meshTriangulated.Cell0DsCoordinates,
									 {});
	
			utilities.ExportSegments("./Cell1Ds.inp",
									 meshTriangulated.Cell0DsCoordinates,
									 meshTriangulated.Cell1DsExtrema,
									 {},
									 {});
	
			utilities.ExportPolygons("./Cell2Ds.inp",
									 meshTriangulated.Cell0DsCoordinates,
									 meshTriangulated.Cell2DsVertices,
									 {},
									 {});
									 
									 //UCDUtilities::ExportPolyhedra
		 }
	}
}