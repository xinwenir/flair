// ********************************************
// Author:  wioletta.kozlowska@cern.ch
// based on the T.Boehlen version from 2011
// Modified wioletta.kozlowska@cern.ch
// Version: 3.1
// Last change: 19/09/2016
// *********************************************
#ifndef ROI_H
#define ROI_H

#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "voxel.h"
#include "vector.h"
#include "eventbin.h"

typedef struct {
	int V[3];
} VECTOR_INT;

typedef struct {
	std::string	Tag;			// Name of the ROI
	int		Id;			// >0 PTV, <0 OAR , 0 non important

	short		RoiId;			// Number of ROI from voxel
	int		VoxelNb;
	double		TargetDose;

	//new
	double		TargetRBE;
	double		TargetLET;
	int		VoxelIndexMax;

	std::vector<int>	VoxelIndex;	// Voxel index in DOSE_GRID
	std::vector<float>	Weight;		// Weight of the voxel-> partially in ROI?

	//only use the following if really needed
	std::vector<int*>	MatrixIndex;	// PB_MATRIX voxel index for each pencil beam
} ROI;

class Roi {
public:
	Vector		Min,Max,Delta, BinNb;		// minimum and maximum point of voxel cube
	double		VoxelVolume;
	int		RoiNb;
	std::string	fileROI;
	std::string	vxlFile;
	std::string	evtFile;
	std::vector<ROI> box;
public:
	Roi() {};
	Roi(std::string filename, float nXmin, float nYmin, float nZmin,
	    std::string evtBin, std::vector<int> numRoi, std::vector<int> isPTV,
	    std::vector<float> TargetDose, std::vector<float> TargetRBE, std::vector<float> TargetLET);
	Roi(std::string filename);
	~Roi();

	int clean(std::vector<ROI> *Struct);
	int readRoi (std::string filename);

	int readVxl (std::string filename, float nXmin, float nYmin, float nZmin, std::string evtBin);

	// TODO correct this
	Vector index2Position(int VoxelIndex);
	int position2Index(Vector Position, int Id);
	int bin2Index(int CurrentBin[3], int Id);
	VECTOR_INT index2Bin(int VoxelIndex);
	bool isInsideRoi(int Bin[3]); //maybe check not only bin also points - see old optimizer
}; // class Roi
//	int manipulateRoi(int RoiBoxId); //we need some manipulation of the target weight
#endif
