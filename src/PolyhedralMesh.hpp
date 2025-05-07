#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolyhedralLibrary {

struct PolyhedralMesh
{
	// VERTICI
    unsigned int NumCell0Ds = 0; // numero di celle 0D
    std::vector<unsigned int> Cell0DsId = {}; // id celle 0D
    Eigen::MatrixXd Cell0DsCoordinates = {}; // coordinate celle 0D
    
    // LATI/SPIGOLI
    unsigned int NumCell1Ds = 0; // numero di celle 1D
    std::vector<unsigned int> Cell1DsId = {}; // id celle 1D
    Eigen::MatrixXi Cell1DsExtrema = {}; // id dei vertici (partenza, arrivo) celle 1D
    
    // FACCE
    unsigned int NumCell2Ds = 0; // numero celle 2D
    std::vector<unsigned int> Cell2DsId = {}; // id celle 2D
    std::vector<std::vector<unsigned int>> Cell2DsVertices = {}; // id dei vertici celle 2D
    std::vector<std::vector<unsigned int>> Cell2DsEdges = {}; // id dei lati celle 2D
    
    // POLIEDRI
    unsigned int NumCell3Ds = 0; // numero celle 3D
    std::vector<unsigned int> Cell3DsId = {}; // id celle 3D
    std::vector<unsigned int> Cell3DsVertices = {}; // id dei vertici celle 3D
    std::vector<unsigned int> Cell3DsEdges = {}; // id dei lati celle 3D
    std::vector<unsigned int> Cell3DsFaces = {}; // id delle facce celle 3D   

};

}
