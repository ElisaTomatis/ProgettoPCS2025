#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <Eigen/Dense>

void ExportParaview(const PolyhedralMesh meshTriangulated){
	
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

        utilities.ExportPolyhedra("./Cell3Ds.inp",
                                 meshTriangulated.Cell0DsCoordinates,
                                 meshTriangulated.Cell3DsVertices,
                                 {},
                                 {});
     }
}
